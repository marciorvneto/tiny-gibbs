#define TINY_LA_IMPLEMENTATION
#define TINY_GIBBS_IMPLEMENTATION
#include "tiny-gibbs.h"
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
int solve(double T, double P, double Pref, double *n, double *nT,
          double *lambda, double *nu, double *total_atoms, int N, int C) {

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

  int result = gibbs_solve_nr(n, nT, lambda, nu, T, P, Pref, gamma, total_atoms,
                              N, C, A, g_db.components, &arena);

  tla_arena_destroy(&arena);
  return result;
}

// WASM entry: solve for adiabatic flame temperature
// Returns T_ad (K), or negative on failure
// Writes equilibrium composition into n[], nT, and gamma_ratio
EMSCRIPTEN_KEEPALIVE
double solve_adiabatic_temperature(double T_feed, double T_guess, double P,
                                   double Pref, double *feed, double *n,
                                   double *nT, double *gamma_ratio,
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

  AdiabaticResult res =
      solve_adiabatic(T_feed, T_guess, P, Pref, feed, total_atoms, N, C, A,
                      g_db.components, &scratch);

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

    gibbs_solve_nr(n, nT, lambda, &nu_val, res.T_ad, P, Pref, gamma,
                   total_atoms, N, C, A, g_db.components, &scratch);

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

  double P = 1;
  double Pref = 1;
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

    int result =
        gibbs_solve_nr(n, &nT, lambda, &nu, T, P, Pref, gamma, total_atoms, N,
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
        solve_adiabatic(T_feed, T_guess, P, Pref, feed_norm, total_atoms, N, C,
                        A, g_db.components, &scratch);

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

      gibbs_solve_nr(n, &nT, lambda, &nu, res.T_ad, P, Pref, gamma, total_atoms,
                     N, C, A, g_db.components, &scratch);

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
