#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tinyla.h>

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

Arena arena_create(size_t capacity);
void arena_destroy(Arena *a);
void *arena_alloc(Arena *a, size_t amount);

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

//================================
//
//   Thermodynamics
//
//================================

#define R_GAS 8.314472

double enthalpy_nondim(DBComponent *sp, double T);
double entropy_nondim(DBComponent *sp, double T);
double gibbs_nondim(DBComponent *sp, double T);

double cp_nondim(DBComponent *sp, double T);

double mixture_enthalpy_over_R(double *n, int N, double T,
                               DBComponent *species);

double mixture_cp_over_R(double *n, int N, double T, DBComponent *species);

// Assuming ideal gas, so that Cp - Cv = R
double mixture_gamma(double *n, double nT, int N, double T,
                     DBComponent *species);

double potential_RT(double n_i, double nT, double T, double P, double Pref,
                    double gamma_i, DBComponent *sp);

void build_jacobian(tla_Matrix *J, double *n, double nT, int N, int C,
                    tla_Matrix *A);

// Build the RHS (negative gradient of Lagrangian)
void build_rhs(tla_Vector *rhs, double *n, double nT, double *lambda, double nu,
               double T, double P, double Pref, double *gamma,
               double *total_atoms, int N, int C, tla_Matrix *A,
               DBComponent *species);

//================================
//
//   Numerical
//
//================================

int gibbs_solve_nr(double *n, double *nT, double *lambda, double *nu, double T,
                   double P, double Pref, double *gamma, double *total_atoms,
                   int N, int C, tla_Matrix *A, DBComponent *species,
                   tla_Arena *scratch);

typedef struct {
  double T_ad;
  int outer_iters;
  int status; // 0 = converged, -1 = inner failed, -2 = outer didn't converge
} AdiabaticResult;

AdiabaticResult solve_adiabatic(double T_feed, double T_guess, double P,
                                double Pref, double *feed, double *total_atoms,
                                int N, int C, tla_Matrix *A,
                                DBComponent *species, tla_Arena *scratch);

#ifdef TINY_GIBBS_IMPLEMENTATION

//================================
//
//   Arena
//
//================================

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

double potential_RT(double n_i, double nT, double T, double P, double Pref,
                    double gamma_i, DBComponent *sp) {
  return gibbs_nondim(sp, T) + log(n_i) - log(nT) + log(gamma_i) +
         log(P / Pref);
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
               double T, double P, double Pref, double *gamma,
               double *total_atoms, int N, int C, tla_Matrix *A,
               DBComponent *species) {

  // dL/dn_i = μ_i/RT + Σ_j λ_j a_ij + ν
  for (int i = 0; i < N; i++) {
    double mu = potential_RT(n[i], nT, T, P, Pref, gamma[i], &species[i]);
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
                   double P, double Pref, double *gamma, double *total_atoms,
                   int N, int C, tla_Matrix *A, DBComponent *species,
                   tla_Arena *scratch) {

  int S = N + 1 + C + 1; // total system size
  double tol = 1e-8;
  int max_iter = 100;

  for (int iter = 0; iter < max_iter; iter++) {
    size_t save = tla_arena_save(scratch);

    // Build RHS and Jacobian
    tla_Vector *rhs = tla_vector_create(scratch, S);
    tla_Matrix *J = tla_matrix_create(scratch, S, S);

    build_rhs(rhs, n, *nT, lambda, *nu, T, P, Pref, gamma, total_atoms, N, C, A,
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

AdiabaticResult solve_adiabatic(double T_feed, double T_guess, double P,
                                double Pref, double *feed, double *total_atoms,
                                int N, int C, tla_Matrix *A,
                                DBComponent *species, tla_Arena *scratch) {

  AdiabaticResult res = {0};
  double tol = 0.5;
  int max_outer = 50;

  double H_react = mixture_enthalpy_over_R(feed, N, T_feed, species);
  double T = T_guess;

  // Working arrays for inner solver
  double n[MAX_COMPONENTS];
  double nT;
  double lambda[5];
  double nu;
  double gamma[MAX_COMPONENTS];

  // 1. INITIALIZE ONCE (Warm Start)
  for (int i = 0; i < N; i++)
    n[i] = 1.0 / N;
  nT = 1.0;
  for (int j = 0; j < C; j++)
    lambda[j] = 1.0;
  nu = 1.0;
  for (int i = 0; i < N; i++)
    gamma[i] = 1.0;

  // Variables for the Secant method
  double T_old = T;
  double f_old = 0.0;

  for (int outer = 0; outer < max_outer; outer++) {

    // Inner solver uses the previous loop's equilibrium as the guess
    int inner_result = gibbs_solve_nr(n, &nT, lambda, &nu, T, P, Pref, gamma,
                                      total_atoms, N, C, A, species, scratch);
    if (inner_result < 0) {
      res.status = -1;
      res.T_ad = T;
      res.outer_iters = outer;
      return res;
    }

    double H_prod = mixture_enthalpy_over_R(n, N, T, species);
    double f = H_prod - H_react;
    double Cp = mixture_cp_over_R(n, N, T, species);

    if (fabs(f / Cp) < tol) {
      res.T_ad = T;
      res.outer_iters = outer + 1;
      res.status = 0;
      return res;
    }

    double dT;
    // 2. SECANT METHOD
    if (outer == 0) {
      // First iteration: Fallback to frozen Cp to generate a second point
      dT = -f / Cp;
    } else {
      double df_dT = (f - f_old) / (T - T_old);

      // Safety net: if slope is too flat or negative (unphysical), use frozen
      // Cp
      if (fabs(df_dT) < 1e-6 || df_dT < 0) {
        dT = -f / Cp;
      } else {
        dT = -f / df_dT;
      }
    }

    // Clamp step to avoid wild jumps
    if (dT > 100.0)
      dT = 100.0;
    if (dT < -100.0)
      dT = -100.0;

    T_old = T;
    f_old = f;

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

#endif
