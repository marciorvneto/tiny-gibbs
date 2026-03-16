#define TINY_LA_IMPLEMENTATION
#include "tinyla.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#define DB_PATH "./db/nasa9_combustion.dat"
#define MAX_LINE_SIZE 1024
#define MAX_ATOMS 128
#define MAX_COMPONENTS 128

//================================
//
//   Arena
//
//================================

typedef struct {
  char *base;
  size_t offset;
  size_t capacity;
} Arena;

Arena arena_create(size_t capacity) {
  Arena a = {0};
  a.capacity = capacity;
  a.base = malloc(capacity);
  return a;
}

void arena_destroy(Arena *a) { free(a->base); }

void *arena_alloc(Arena *a, size_t amount) {
  size_t aligned = (amount + 7) & ~7;
  if (a->offset + aligned > a->capacity) {
    return NULL;
  }
  void *addr = a->base + a->offset;
  a->offset += aligned;
  return addr;
}

//================================
//
//   IO
//
//================================

const char *read_db_file(Arena *a, const char *path) {
  FILE *db_file = fopen(path, "r");
  if (db_file == NULL) {
    printf("Could not open the database file\n");
    exit(1);
  }

  if (fseek(db_file, 0, SEEK_END) != 0) {
    fclose(db_file);
    return NULL;
  }

  long size = ftell(db_file);

  rewind(db_file);

  char *buffer = arena_alloc(a, (size_t)size + 1); // Null-terminated
  if (buffer == NULL) {
    printf("Arena ran out of memory when reading database.\n");
    exit(1);
  }

  size_t bytes_read = fread(buffer, 1, size, db_file);
  buffer[bytes_read] = '\0';

  fclose(db_file);

  return buffer;
}

//================================
//
//   DB
//
//================================

typedef struct {
  char symbol[4]; // "C", "H", "O", "N", "Ar"
  int count;
} AtomEntry;

typedef struct {
  char name[16];
  int n_elements;
  double mol_weight;
  double hf298; // J/mol
  double T_low, T_mid, T_high;
  AtomEntry atoms[5]; // max 5 elements (C,H,O,N,Ar)
  double low[9];      // NASA-9 coeffs, low-T range
  double high[9];     // NASA-9 coeffs, high-T range
} DBComponent;

typedef struct {
  DBComponent components[MAX_COMPONENTS];
  size_t num_components;
} DBData;

//================================
//
//   DB Parser
//
//================================

typedef struct {
  const char *string;
  size_t pointer;
  size_t line;
} Parser;

typedef struct {
  const char *start;
  size_t length;
  size_t line;
} LineView;

Parser create_parser(const char *string) {
  Parser p = {0};
  p.line = 1;
  p.string = string;
  return p;
}

// Returns the number of bytes read
int read_line(Parser *p, LineView *line) {
  if (p->string[p->pointer] == '\0')
    return 0;

  size_t start = p->pointer;
  size_t end = start;

  while (p->string[end] != '\n' && p->string[end] != '\0') {
    end++;
  }

  line->start = p->string + start;
  line->length = end - start;
  line->line = p->line;

  p->pointer = end; // Default
  if (p->string[end] == '\n') {
    p->line++;
    p->pointer = end + 1;
  }
  if (p->string[end] == '\0') {
    return 0;
  }

  return end - start;
}

void parse_component_atoms(const char *line_str, DBComponent *c) {

  const char *cursor = line_str;
  int n;
  for (size_t i = 0; i < c->n_elements; i++) {
    sscanf(cursor, "%3s %d%n", c->atoms[i].symbol, &c->atoms[i].count, &n);
    cursor += n;
  }
}

void parse_db(Parser *p, DBData *db) {
  LineView line = {0};
  char line_str[MAX_LINE_SIZE];
  int num_components = -1;

  // Skip comments and get to the number of components
  while (p->string[p->pointer] != '\0') {
    read_line(p, &line);
    if (line.start[0] == '#')
      continue;
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    num_components = atoi(line_str);
    break;
  }

  for (int i = 0; i < num_components; i++) {
    DBComponent *c = &db->components[db->num_components];

    // Line 1: name  n_elements  mol_weight  hf298  T_low  T_mid  T_high
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%15s %d %lf %lf %lf %lf %lf", c->name, &c->n_elements,
           &c->mol_weight, &c->hf298, &c->T_low, &c->T_mid, &c->T_high);

    // Line 2: element pairs
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    parse_component_atoms(line_str, c);

    // Line 3: low-T coefficients
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &c->low[0],
           &c->low[1], &c->low[2], &c->low[3], &c->low[4], &c->low[5],
           &c->low[6], &c->low[7], &c->low[8]);

    // Line 4: high-T coefficients
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &c->high[0],
           &c->high[1], &c->high[2], &c->high[3], &c->high[4], &c->high[5],
           &c->high[6], &c->high[7], &c->high[8]);

    db->num_components++;
  }
}

//================================
//
//   Thermodynamics
//
//================================

#define R_GAS 8.314472

double enthalpy_nondim(DBComponent *sp, double T) {
  double *c = (T <= sp->T_mid) ? sp->low : sp->high;
  double T2 = T * T, T3 = T2 * T, T4 = T3 * T, lnT = log(T);

  return -c[0] / (T2) + c[1] * lnT / T + c[2] + c[3] * T / 2.0 +
         c[4] * T2 / 3.0 + c[5] * T3 / 4.0 + c[6] * T4 / 5.0 + c[7] / T;
}

double entropy_nondim(DBComponent *sp, double T) {
  double *c = (T <= sp->T_mid) ? sp->low : sp->high;
  double T2 = T * T, T3 = T2 * T, T4 = T3 * T, lnT = log(T);

  return -c[0] / (2.0 * T2) - c[1] / T + c[2] * lnT + c[3] * T +
         c[4] * T2 / 2.0 + c[5] * T3 / 3.0 + c[6] * T4 / 4.0 + c[8];
}

double gibbs_nondim(DBComponent *sp, double T) {
  double *c = (T <= sp->T_mid) ? sp->low : sp->high;
  double T2 = T * T, T3 = T2 * T, T4 = T3 * T, lnT = log(T);

  double H_RT = -c[0] / (T2) + c[1] * lnT / T + c[2] + c[3] * T / 2.0 +
                c[4] * T2 / 3.0 + c[5] * T3 / 4.0 + c[6] * T4 / 5.0 + c[7] / T;

  double S_R = -c[0] / (2.0 * T2) - c[1] / T + c[2] * lnT + c[3] * T +
               c[4] * T2 / 2.0 + c[5] * T3 / 3.0 + c[6] * T4 / 4.0 + c[8];

  return H_RT - S_R;
}

double cp_nondim(DBComponent *sp, double T) {
  double *c = (T <= sp->T_mid) ? sp->low : sp->high;
  double T2 = T * T, T3 = T2 * T, T4 = T3 * T;

  return c[0] / (T2) + c[1] / T + c[2] + c[3] * T + c[4] * T2 + c[5] * T3 +
         c[6] * T4;
}

double mixture_enthalpy_over_R(double *n, int N, double T,
                               DBComponent *species) {
  double H = 0.0;
  for (int i = 0; i < N; i++) {
    H += n[i] * enthalpy_nondim(&species[i], T);
  }
  return H * T; // H_mix / R  [K·mol]
}

double mixture_cp_over_R(double *n, int N, double T, DBComponent *species) {
  double Cp = 0.0;
  for (int i = 0; i < N; i++) {
    Cp += n[i] * cp_nondim(&species[i], T);
  }
  return Cp; // Cp_mix / R
}

// Assuming ideal gas, so that Cp - Cv = R
double mixture_gamma(double *n, double nT, int N, double T,
                     DBComponent *species) {
  double Cp_over_R = mixture_cp_over_R(n, N, T, species);
  double Cv_over_R = Cp_over_R - nT;
  return Cp_over_R / Cv_over_R;
}

double potential_RT(double n_i, double nT, double T, double gamma_i,
                    DBComponent *sp) {
  return gibbs_nondim(sp, T) + log(n_i) - log(nT) + log(gamma_i);
}

// TODO: Include the derivation in the repo
// Build the KKT Jacobian: size (N+1+C+1) x (N+1+C+1)
//  [ diag(1/n)   -1/nT    A      1  ]
//  [ -1/nT       Σn/nT²   0     -1  ]
//  [  A^T         0        0      0  ]
//  [  1^T        -1        0      0  ]
// where A is the atom matrix (N x C)
void build_jacobian(tla_Matrix *J, double *n, double nT, int N, int C,
                    tla_Matrix *A) {
  int S = N + 1 + C + 1;

  // Zero everything
  for (int i = 0; i < S; i++)
    for (int j = 0; j < S; j++)
      tla_matrix_set_value(J, i, j, 0.0);

  // Block: diag(1/n_i)
  for (int i = 0; i < N; i++)
    tla_matrix_set_value(J, i, i, 1.0 / n[i]);

  // Block: -1/nT column and row
  for (int i = 0; i < N; i++) {
    tla_matrix_set_value(J, i, N, -1.0 / nT); // col N
    tla_matrix_set_value(J, N, i, -1.0 / nT); // row N
  }

  // Block: nT diagonal
  double sum_n = 0;
  for (int i = 0; i < N; i++)
    sum_n += n[i];
  tla_matrix_set_value(J, N, N, sum_n / (nT * nT));

  // Block: A (atom matrix) and A^T
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < C; j++) {
      double a_ij = tla_matrix_get_value(A, i, j);
      tla_matrix_set_value(J, i, N + 1 + j, a_ij); // A
      tla_matrix_set_value(J, N + 1 + j, i, a_ij); // A^T
    }
  }

  // Block: ν column and row (ones)
  int nu_col = N + 1 + C;
  for (int i = 0; i < N; i++) {
    tla_matrix_set_value(J, i, nu_col, 1.0);
    tla_matrix_set_value(J, nu_col, i, 1.0);
  }

  // nT-ν entries
  tla_matrix_set_value(J, N, nu_col, -1.0);
  tla_matrix_set_value(J, nu_col, N, -1.0);
}

// Build the RHS (negative gradient of Lagrangian)
void build_rhs(tla_Vector *rhs, double *n, double nT, double *lambda, double nu,
               double T, double *gamma, double *total_atoms, int N, int C,
               tla_Matrix *A, DBComponent *species) {

  // dL/dn_i = μ_i/RT + Σ_j λ_j a_ij + ν
  for (int i = 0; i < N; i++) {
    double mu = potential_RT(n[i], nT, T, gamma[i], &species[i]);
    double atom_term = 0;
    for (int j = 0; j < C; j++)
      atom_term += lambda[j] * tla_matrix_get_value(A, i, j);
    tla_vector_set_value(rhs, i, -(mu + atom_term + nu));
  }

  // dL/dnT = -Σn/nT - ν
  double sum_n = 0;
  for (int i = 0; i < N; i++)
    sum_n += n[i];
  tla_vector_set_value(rhs, N, -(-sum_n / nT - nu));

  // dL/dλ_j = Σ_i n_i a_ij - b_j
  for (int j = 0; j < C; j++) {
    double balance = -total_atoms[j];
    for (int i = 0; i < N; i++)
      balance += n[i] * tla_matrix_get_value(A, i, j);
    tla_vector_set_value(rhs, N + 1 + j, -balance);
  }

  // dL/dν = Σn - nT
  tla_vector_set_value(rhs, N + 1 + C, -(sum_n - nT));
}

//================================
//
//   Numerical
//
//================================

int gibbs_solve_nr(double *n, double *nT, double *lambda, double *nu, double T,
                   double *gamma, double *total_atoms, int N, int C,
                   tla_Matrix *A, DBComponent *species, tla_Arena *scratch) {

  int S = N + 1 + C + 1; // total system size
  double tol = 1e-8;
  int max_iter = 100;

  for (int iter = 0; iter < max_iter; iter++) {
    size_t save = tla_arena_save(scratch);

    // Build RHS and Jacobian
    tla_Vector *rhs = tla_vector_create(scratch, S);
    tla_Matrix *J = tla_matrix_create(scratch, S, S);

    build_rhs(rhs, n, *nT, lambda, *nu, T, gamma, total_atoms, N, C, A,
              species);
    build_jacobian(J, n, *nT, N, C, A);

    // Check convergence (infinity norm of rhs, which is -F)
    double err = 0;
    for (int i = 0; i < S; i++) {
      double v = fabs(tla_vector_get_value(rhs, i));
      if (v > err)
        err = v;
    }
    if (err < tol) {
      tla_arena_restore(scratch, save);
      return iter;
    }

    // Solve J * dx = rhs  (rhs is already -F)
    tla_Matrix *aug = tla_matrix_append_column(scratch, J, rhs);
    int code;
    tla_Vector *dx = gauss_solve_new(scratch, aug, &code);
    if (code != 0) {
      tla_arena_restore(scratch, save);
      return -1; // singular
    }

    // Unpack
    double *dn = dx->values;            // [0..N-1]
    double dnT = dx->values[N];         // [N]
    double *dlamb = dx->values + N + 1; // [N+1..N+C]
    double dnu = dx->values[N + 1 + C]; // [N+1+C]

    // IMPORTANT: Do not skip this step!
    // Damping: fraction to boundary
    double alpha = 1.0;
    double tau = 0.99;
    for (int i = 0; i < N; i++) {
      if (dn[i] < 0.0) {
        double max_step = -n[i] / dn[i];
        if (tau * max_step < alpha)
          alpha = tau * max_step;
      }
    }
    if (dnT < 0.0) {
      double max_step = -(*nT) / dnT;
      if (tau * max_step < alpha)
        alpha = tau * max_step;
    }

    // Apply step
    for (int i = 0; i < N; i++)
      n[i] += alpha * dn[i];
    *nT += alpha * dnT;
    for (int j = 0; j < C; j++)
      lambda[j] += alpha * dlamb[j];
    *nu += alpha * dnu;

    tla_arena_restore(scratch, save);
  }
  return -2; // didn't converge
}

typedef struct {
  double T_ad;
  int outer_iters;
  int status; // 0 = converged, -1 = inner failed, -2 = outer didn't converge
} AdiabaticResult;

AdiabaticResult solve_adiabatic(double T_feed, double T_guess, double *feed,
                                double *total_atoms, int N, int C,
                                tla_Matrix *A, DBComponent *species,
                                tla_Arena *scratch) {

  AdiabaticResult res = {0};
  double tol = 0.5; // 0.5 K tolerance on T
  int max_outer = 50;

  // Compute reactant enthalpy at feed temperature
  double H_react = mixture_enthalpy_over_R(feed, N, T_feed, species);

  double T = T_guess;

  // Working arrays for inner solver — re-initialized each outer iteration
  double n[MAX_COMPONENTS];
  double nT;
  double lambda[5];
  double nu;
  double gamma[MAX_COMPONENTS];

  for (int outer = 0; outer < max_outer; outer++) {
    // Re-initialize inner solver state for this T
    for (int i = 0; i < N; i++)
      n[i] = 1.0 / N;
    nT = 1.0;
    for (int j = 0; j < C; j++)
      lambda[j] = 1.0;
    nu = 1.0;
    for (int i = 0; i < N; i++)
      gamma[i] = 1.0;

    // Solve equilibrium composition at current T
    int inner_result = gibbs_solve_nr(n, &nT, lambda, &nu, T, gamma,
                                      total_atoms, N, C, A, species, scratch);
    if (inner_result < 0) {
      res.status = -1;
      res.T_ad = T;
      res.outer_iters = outer;
      return res;
    }

    // Evaluate energy residual: f(T) = H_prod(T) - H_react(T_feed)
    double H_prod = mixture_enthalpy_over_R(n, N, T, species);
    double f = H_prod - H_react;

    // Check convergence
    // |f| < tol * Cp gives us sub-Kelvin accuracy
    double Cp = mixture_cp_over_R(n, N, T, species);
    if (fabs(f / Cp) < tol) {
      res.T_ad = T;
      res.outer_iters = outer + 1;
      res.status = 0;
      return res;
    }

    // Newton step: dT = -f / Cp
    double dT = -f / Cp;

    // Clamp step to avoid wild jumps
    if (dT > 500.0)
      dT = 500.0;
    if (dT < -500.0)
      dT = -500.0;

    // Keep T in a physically reasonable range
    T += dT;
    if (T < 300.0)
      T = 300.0;
    if (T > 6000.0)
      T = 6000.0;
  }

  res.T_ad = T;
  res.outer_iters = max_outer;
  res.status = -2;
  return res;
}

//================================
//
//   WASM API
//
//================================

static DBData g_db;
static const char *g_elements[] = {"C", "H", "O", "N", "Ar"};

EMSCRIPTEN_KEEPALIVE
void init() {
  Arena a = arena_create(1024 * 1024);
  const char *contents = read_db_file(&a, DB_PATH);
  Parser p = create_parser(contents);
  parse_db(&p, &g_db);
  arena_destroy(&a);
}

EMSCRIPTEN_KEEPALIVE
int solve(double T, double *n, double *nT, double *lambda, double *nu,
          double *total_atoms, int N, int C) {

  tla_Arena arena = tla_arena_create(4 * 1024 * 1024);

  tla_Matrix *A = tla_matrix_of_value(&arena, N, C, 0.0);
  for (int i = 0; i < N; i++) {
    DBComponent *sp = &g_db.components[i];
    for (int k = 0; k < sp->n_elements; k++) {
      for (int j = 0; j < C; j++) {
        if (strcmp(sp->atoms[k].symbol, g_elements[j]) == 0) {
          tla_matrix_set_value(A, i, j, (double)sp->atoms[k].count);
          break;
        }
      }
    }
  }

  double *gamma = tla_arena_alloc(&arena, N * sizeof(double));
  for (int i = 0; i < N; i++)
    gamma[i] = 1.0;

  int result = gibbs_solve_nr(n, nT, lambda, nu, T, gamma, total_atoms, N, C, A,
                              g_db.components, &arena);

  tla_arena_destroy(&arena);
  return result;
}

// WASM entry: solve for adiabatic flame temperature
// Returns T_ad (K), or negative on failure
// Writes equilibrium composition into n[], nT, and gamma_ratio
EMSCRIPTEN_KEEPALIVE
double solve_adiabatic_temperature(double T_feed, double T_guess, double *feed,
                                   double *n, double *nT, double *gamma_ratio,
                                   double *total_atoms, int N, int C,
                                   int *iterations) {

  tla_Arena arena = tla_arena_create(4 * 1024 * 1024);

  // Build atom matrix
  tla_Matrix *A = tla_matrix_of_value(&arena, N, C, 0.0);
  for (int i = 0; i < N; i++) {
    DBComponent *sp = &g_db.components[i];
    for (int k = 0; k < sp->n_elements; k++) {
      for (int j = 0; j < C; j++) {
        if (strcmp(sp->atoms[k].symbol, g_elements[j]) == 0) {
          tla_matrix_set_value(A, i, j, (double)sp->atoms[k].count);
          break;
        }
      }
    }
  }

  tla_Arena scratch = tla_arena_create(4 * 1024 * 1024);

  AdiabaticResult res = solve_adiabatic(T_feed, T_guess, feed, total_atoms, N,
                                        C, A, g_db.components, &scratch);

  // If converged, do one final equilibrium solve at T_ad to populate n[] and nT
  if (res.status == 0) {
    double lambda[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double nu_val = 1.0;
    double gamma[MAX_COMPONENTS];

    for (int i = 0; i < N; i++)
      n[i] = 1.0 / N;
    *nT = 1.0;
    for (int i = 0; i < N; i++)
      gamma[i] = 1.0;

    gibbs_solve_nr(n, nT, lambda, &nu_val, res.T_ad, gamma, total_atoms, N, C,
                   A, g_db.components, &scratch);

    // Compute heat capacity ratio at equilibrium
    *gamma_ratio = mixture_gamma(n, *nT, N, res.T_ad, g_db.components);
    *iterations = res.outer_iters;
  }

  tla_arena_destroy(&scratch);
  tla_arena_destroy(&arena);

  return (res.status == 0) ? res.T_ad : -(double)(res.status);
}

int main() {
  init();

  int N = g_db.num_components;
  int C = 5;

  tla_Arena la = tla_arena_create(1024 * 1024);
  tla_Matrix *A = tla_matrix_of_value(&la, N, C, 0.0);
  for (int i = 0; i < N; i++) {
    DBComponent *sp = &g_db.components[i];
    for (int k = 0; k < sp->n_elements; k++) {
      for (int j = 0; j < C; j++) {
        if (strcmp(sp->atoms[k].symbol, g_elements[j]) == 0) {
          tla_matrix_set_value(A, i, j, (double)sp->atoms[k].count);
          break;
        }
      }
    }
  }

  // CH4/air stoichiometric
  double feed[MAX_COMPONENTS] = {0};
  int idx_ch4 = -1, idx_o2 = -1, idx_n2 = -1, idx_ar = -1;
  for (int i = 0; i < N; i++) {
    if (strcmp(g_db.components[i].name, "CH4") == 0)
      idx_ch4 = i;
    if (strcmp(g_db.components[i].name, "O2") == 0)
      idx_o2 = i;
    if (strcmp(g_db.components[i].name, "N2") == 0)
      idx_n2 = i;
    if (strcmp(g_db.components[i].name, "Ar") == 0)
      idx_ar = i;
  }

  feed[idx_ch4] = 1.0;
  feed[idx_o2] = 2.0;
  feed[idx_n2] = 7.52;
  feed[idx_ar] = 0.09;

  // Compute atom balances (normalized per mole of mixture)
  double total_atoms[5] = {0};
  double total_feed = 0;
  for (int i = 0; i < N; i++)
    total_feed += feed[i];
  for (int j = 0; j < C; j++)
    for (int i = 0; i < N; i++)
      total_atoms[j] += feed[i] * tla_matrix_get_value(A, i, j);
  for (int j = 0; j < C; j++)
    total_atoms[j] /= total_feed;

  // Normalize feed too (per mole of mixture)
  double feed_norm[MAX_COMPONENTS] = {0};
  for (int i = 0; i < N; i++)
    feed_norm[i] = feed[i] / total_feed;

  //===================================
  // Test 1: Fixed-T equilibrium at 2000 K
  //===================================
  printf("=== Fixed-T Equilibrium (T = 2000 K) ===\n\n");
  {
    double n[MAX_COMPONENTS];
    for (int i = 0; i < N; i++)
      n[i] = 1.0 / N;
    double nT = 1.0;
    double lambda[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double nu = 1.0;
    double gamma[MAX_COMPONENTS];
    for (int i = 0; i < N; i++)
      gamma[i] = 1.0;

    double T = 2000.0;
    tla_Arena scratch = tla_arena_create(4 * 1024 * 1024);

    int result = gibbs_solve_nr(n, &nT, lambda, &nu, T, gamma, total_atoms, N,
                                C, A, g_db.components, &scratch);

    if (result >= 0) {
      printf("Converged in %d iterations at T = %.1f K\n\n", result, T);
      printf("%-8s %12s %12s\n", "Species", "n_i", "x_i");
      printf("--------------------------------------\n");
      for (int i = 0; i < N; i++) {
        if (n[i] > 1e-6)
          printf("%-8s %12.6f %12.6f\n", g_db.components[i].name, n[i],
                 n[i] / nT);
      }
      printf("\nTotal moles: %.6f\n", nT);

      double gam = mixture_gamma(n, nT, N, T, g_db.components);
      printf("Heat capacity ratio (gamma): %.4f\n", gam);
    } else {
      printf("Failed (code %d)\n", result);
    }

    tla_arena_destroy(&scratch);
  }

  //===================================
  // Test 2: Adiabatic flame temperature
  //===================================
  printf("\n=== Adiabatic Flame Temperature ===\n\n");
  {
    double T_feed = 298.15;
    double T_guess = 2000.0;

    tla_Arena scratch = tla_arena_create(4 * 1024 * 1024);

    AdiabaticResult res =
        solve_adiabatic(T_feed, T_guess, feed_norm, total_atoms, N, C, A,
                        g_db.components, &scratch);

    if (res.status == 0) {
      printf("Converged in %d outer iterations\n", res.outer_iters);
      printf("Adiabatic flame temperature: %.1f K\n\n", res.T_ad);

      // Final equilibrium at T_ad for display
      double n[MAX_COMPONENTS];
      for (int i = 0; i < N; i++)
        n[i] = 1.0 / N;
      double nT = 1.0;
      double lambda[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
      double nu = 1.0;
      double gamma[MAX_COMPONENTS];
      for (int i = 0; i < N; i++)
        gamma[i] = 1.0;

      gibbs_solve_nr(n, &nT, lambda, &nu, res.T_ad, gamma, total_atoms, N, C, A,
                     g_db.components, &scratch);

      printf("%-8s %12s %12s\n", "Species", "n_i", "x_i");
      printf("--------------------------------------\n");
      for (int i = 0; i < N; i++) {
        if (n[i] > 1e-6)
          printf("%-8s %12.6f %12.6f\n", g_db.components[i].name, n[i],
                 n[i] / nT);
      }
      printf("\nTotal moles: %.6f\n", nT);

      double gam = mixture_gamma(n, nT, N, res.T_ad, g_db.components);
      printf("Heat capacity ratio (gamma): %.4f\n", gam);
    } else if (res.status == -1) {
      printf("Failed: inner solver failed at T = %.1f K (iter %d)\n", res.T_ad,
             res.outer_iters);
    } else {
      printf(
          "Failed: outer loop did not converge (T = %.1f K after %d iters)\n",
          res.T_ad, res.outer_iters);
    }

    tla_arena_destroy(&scratch);
  }

  tla_arena_destroy(&la);
  return 0;
}
