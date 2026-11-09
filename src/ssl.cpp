/*******************************************************************************
 * SVMlin
 * Copyright (c) 2006 Vikas Sindhwani at the University of Chicago.
 * Adapted to Percolator by Lukas Käll at the University of Washington
 * Sped up by John Halloran at the University of California, Davis, as detailed in:
 ******************************
 * A Matter of Time: Faster Percolator Analysis via Efficient SVM Learning for 
 * Large-Scale Proteomics
 * John T. Halloran and David M. Rocke
 * Journal of Proteome Research 2018 17 (5), 1978-1982
 *******************************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <ctype.h>
using namespace std;
#include "Globals.h"
#include "ssl.h"

#include <stdarg.h>
#include <cstring>
#include "Timer.h"

extern "C" {
  extern double ddot_(int *, double *, int *, double *, int *); // compute the dot product of two vectors
  extern int daxpy_(int *, double *, double *, int *, double *, int *); // compute y := alpha * x + y
  extern int dscal_(int *, double *, double *, int *); // Compute y := alpha * y
}

#define VERBOSE 1
#define LOG2(x) 1.4426950408889634*log(x)
// for compatibility issues, not using log2

/* Buffers of one L2_SVM_MFN solve, sized for all m examples and reused by
   all its CGLS calls and line searches. */
SolverWorkspace::SolverWorkspace(int m, int n)
    : z(static_cast<size_t>(m)), q(static_cast<size_t>(m)),
      r(static_cast<size_t>(n)), p(static_cast<size_t>(n)),
      deltas(static_cast<size_t>(m)),
      gram(static_cast<size_t>(n) * static_cast<size_t>(n)) {}

/* q[i] = x_{J[i]}' p for the active rows (length n, intercept included).
   Gives the same result as the reference BLAS dgemv with trans='T', alpha=1
   and beta=0 on a copy of the active rows, which this replaces: every row is
   summed in the same order, but four rows are processed together so that
   their additions overlap instead of waiting on each other. */
static void activeRowsTimesVector(int active, const int* J, const AlgIn& data,
                                  const double* p, double* q) {
  const int n = data.n;
  int i = 0;
  for (; i + 4 <= active; i += 4) {
    const double* x0 = data.row(J[i]);
    const double* x1 = data.row(J[i + 1]);
    const double* x2 = data.row(J[i + 2]);
    const double* x3 = data.row(J[i + 3]);
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
    for (int k = 0; k < n; k++) {
      t0 += x0[k] * p[k];
      t1 += x1[k] * p[k];
      t2 += x2[k] * p[k];
      t3 += x3[k] * p[k];
    }
    // dgemv adds the dot product to the zeroed output vector
    q[i] = 0.0 + t0;
    q[i + 1] = 0.0 + t1;
    q[i + 2] = 0.0 + t2;
    q[i + 3] = 0.0 + t3;
  }
  for (; i < active; i++) {
    const double* x = data.row(J[i]);
    double t = 0.0;
    for (int k = 0; k < n; k++) {
      t += x[k] * p[k];
    }
    q[i] = 0.0 + t;
  }
}

double cglsFun1(int active, const int* J, const double* Y,
                const AlgIn& data, double* q,
                double* p, double cpos, double cneg){
  double omega_q = 0.0;
  int i = 0;

  activeRowsTimesVector(active, J, data, p, q);

  for (i = 0; i < active; i++) {
    omega_q += ((Y[J[i]]==1)? cpos : cneg) * (q[i]) * (q[i]);
  }

  return(omega_q);
}

/* r += z[i] * x_i for one row, like the reference BLAS daxpy, which skips
   the row when z[i] is zero. */
static inline void addScaledRow(double zi, const double* x, int n, double* r) {
  if (zi != 0.0) {
    for (int k = 0; k < n; k++) {
      r[k] += zi * x[k];
    }
  }
}

void cglsFun2(int active, const int* J, const double* Y,
              const AlgIn& data, double* q,
              double* o, double* z, double* r, 
              double cpos, double cneg){
  const int n = data.n;
  int i = 0;
  // Rows are added to r four at a time, in the same order per element as
  // one daxpy per row, so that r is loaded and stored once per four rows.
  for (; i + 4 <= active; i += 4) {
    for (int j = i; j < i + 4; j++) {
      o[J[j]] += q[j];
      z[j] -= ((Y[J[j]]==1)? cpos : cneg) * q[j];
    }
    const double* x0 = data.row(J[i]);
    const double* x1 = data.row(J[i + 1]);
    const double* x2 = data.row(J[i + 2]);
    const double* x3 = data.row(J[i + 3]);
    const double z0 = z[i], z1 = z[i + 1], z2 = z[i + 2], z3 = z[i + 3];
    if (z0 != 0.0 && z1 != 0.0 && z2 != 0.0 && z3 != 0.0) {
      for (int k = 0; k < n; k++) {
        double rk = r[k];
        rk += z0 * x0[k];
        rk += z1 * x1[k];
        rk += z2 * x2[k];
        rk += z3 * x3[k];
        r[k] = rk;
      }
    } else {
      addScaledRow(z0, x0, n, r);
      addScaledRow(z1, x1, n, r);
      addScaledRow(z2, x2, n, r);
      addScaledRow(z3, x3, n, r);
    }
  }
  for (; i < active; i++) {
    o[J[i]] += q[i];
    z[i] -= ((Y[J[i]]==1)? cpos : cneg) * q[i];
    addScaledRow(z[i], data.row(J[i]), n, r);
  }
}

/* Dot product with the same summation order as the reference BLAS ddot:
   blocks of five products are summed before they are added. */
#define DDOT_BLOCK(x, y, i)                                         \
  ((x)[i] * (y)[i] + (x)[(i) + 1] * (y)[(i) + 1] +                  \
   (x)[(i) + 2] * (y)[(i) + 2] + (x)[(i) + 3] * (y)[(i) + 3] +      \
   (x)[(i) + 4] * (y)[(i) + 4])

/* o[rows[i]] = x_{rows[i]}' w + w[n0] for the given rows, with the n0
   features of each row and the intercept w[n0]. Same result as one ddot per
   row over the features, but four rows are processed together so that their
   additions and memory loads overlap. */
static void scoreRows(const int* rows, int numRows, const AlgIn& data, int n0,
                      const double* w, double* o) {
  int i = 0;
  for (; i + 4 <= numRows; i += 4) {
    const double* x0 = data.row(rows[i]);
    const double* x1 = data.row(rows[i + 1]);
    const double* x2 = data.row(rows[i + 2]);
    const double* x3 = data.row(rows[i + 3]);
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
    int k = 0;
    for (; k < n0 - 4; k += 5) {
      t0 += DDOT_BLOCK(x0, w, k);
      t1 += DDOT_BLOCK(x1, w, k);
      t2 += DDOT_BLOCK(x2, w, k);
      t3 += DDOT_BLOCK(x3, w, k);
    }
    for (; k < n0; k++) {
      t0 += x0[k] * w[k];
      t1 += x1[k] * w[k];
      t2 += x2[k] * w[k];
      t3 += x3[k] * w[k];
    }
    o[rows[i]] = t0 + w[n0];
    o[rows[i + 1]] = t1 + w[n0];
    o[rows[i + 2]] = t2 + w[n0];
    o[rows[i + 3]] = t3 + w[n0];
  }
  for (; i < numRows; i++) {
    const double* x = data.row(rows[i]);
    double t = 0.0;
    int k = 0;
    for (; k < n0 - 4; k += 5) {
      t += DDOT_BLOCK(x, w, k);
    }
    for (; k < n0; k++) {
      t += x[k] * w[k];
    }
    o[rows[i]] = t + w[n0];
  }
}
#undef DDOT_BLOCK

/* Solves the same least squares problem as CGLS exactly, by Cholesky
   factorisation of the normal equations
     (lambda I + sum_i C_i x_i x_i') beta = sum_i C_i Y_i x_i
   over the active examples. The n x n system (n = number of features plus
   one) is built in a single pass over the active rows, whereas CGLS needs
   two passes per iteration. Updates the outputs of the active examples. */
int directLeastSquares(const AlgIn& data, const double lambda, const int* J,
                       int active, double* beta, double* o, double cpos,
                       double cneg, SolverWorkspace& workspace) {
  const int n = data.n;
  const double* Y = data.Y.data();
  double* G = workspace.gram.data();  // upper triangle, row-major
  double* b = workspace.r.data();
  std::fill(G, G + static_cast<size_t>(n) * n, 0.0);
  std::fill(b, b + n, 0.0);
  for (int i = 0; i < active; i++) {
    const int ii = J[i];
    const double* x = data.row(ii);
    const double c = (Y[ii] == 1) ? cpos : cneg;
    for (int j = 0; j < n; j++) {
      const double cx = c * x[j];
      b[j] += cx * Y[ii];
      double* Gj = G + static_cast<size_t>(j) * n;
      for (int k = j; k < n; k++) {
        Gj[k] += cx * x[k];
      }
    }
  }
  // Cholesky factorisation G = L L', with L stored in the lower triangle
  for (int j = 0; j < n; j++) {
    G[static_cast<size_t>(j) * n + j] += lambda;
    for (int k = j + 1; k < n; k++) {
      G[static_cast<size_t>(k) * n + j] = G[static_cast<size_t>(j) * n + k];
    }
  }
  for (int j = 0; j < n; j++) {
    double* Lj = G + static_cast<size_t>(j) * n;
    double s = Lj[j];
    for (int k = 0; k < j; k++) s -= Lj[k] * Lj[k];
    if (!(s > 0.0)) return 0;  // not positive definite
    const double ljj = std::sqrt(s);
    Lj[j] = ljj;
    for (int i = j + 1; i < n; i++) {
      double* Li = G + static_cast<size_t>(i) * n;
      double t = Li[j];
      for (int k = 0; k < j; k++) t -= Li[k] * Lj[k];
      Li[j] = t / ljj;
    }
  }
  // forward and back substitution
  for (int i = 0; i < n; i++) {
    const double* Li = G + static_cast<size_t>(i) * n;
    double t = b[i];
    for (int k = 0; k < i; k++) t -= Li[k] * beta[k];
    beta[i] = t / Li[i];
  }
  for (int i = n - 1; i >= 0; i--) {
    double t = beta[i];
    for (int k = i + 1; k < n; k++) t -= G[static_cast<size_t>(k) * n + i] * beta[k];
    beta[i] = t / G[static_cast<size_t>(i) * n + i];
  }
  for (int i = 0; i < active; i++) {
    const int ii = J[i];
    const double* x = data.row(ii);
    double t = 0.0;
    for (int k = 0; k < n; k++) t += x[k] * beta[k];
    o[ii] = t;
  }
  return 1;
}

int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
         const double epsilon, const int* J, int active,
         double* beta, double* o,
         double cpos, double cneg, SolverWorkspace& workspace) {
  if (VERBOSE_CGLS) {
    cout << "CGLS starting..." << endl;
  }
  /* Disassemble the structures */
  Timer tictoc;
  const double* Y = data.Y.data();
  int n = data.n;
  // initialize z
  double* z = workspace.z.data();
  double* q = workspace.q.data();
  int ii = 0;
  int i;
  int n0 = n-1;
  int inc = 1;
  double one = 1;
  double negLambda = -lambda;
  double* r = workspace.r.data();
  for (i = n; i--;) {
    r[i] = 0.0;
  }
  for (i = 0; i < active; i++) {
    ii = J[i];
    z[i] = ((Y[ii]==1)? cpos : cneg) * (Y[ii] - o[ii]);
    addScaledRow(z[i], data.row(ii), n, r);
  }
  double* p = workspace.p.data();
  daxpy_(&n, &negLambda, beta, &inc, r, &inc);
  memcpy(p, r, sizeof(double)*static_cast<std::size_t>(n));
  double omega1 = ddot_(&n, r, &inc, r, &inc);
  double omega_p = omega1;
  double omega_q = 0.0;
  double inv_omega2 = 1 / omega1;
  double scale = 0.0;
  double omega_z = 0.0;
  double gamma = 0.0;
  int cgiter = 0;
  int optimality = 0;
  double epsilon2 = epsilon * epsilon;
  // iterate
  while (cgiter < cgitermax) {
    cgiter++;
    omega_q = cglsFun1(active, J, Y, data, q, p, cpos, cneg);
    gamma = omega1 / (lambda * omega_p + omega_q);
    inv_omega2 = 1 / omega1;

    memcpy(r, beta, sizeof(double)*static_cast<std::size_t>(n));
    dscal_(&n, &negLambda, r, &inc);

    daxpy_(&n, &gamma, p, &inc, beta, &inc);
    dscal_(&active, &gamma, q, &inc);

    cglsFun2(active, J, Y, data, q, o, z, r, cpos, cneg);

    omega_z = ddot_(&active, z, &inc, z, &inc);
    omega1 = ddot_(&n, r, &inc, r, &inc);

    if (VERBOSE_CGLS) {
      cout << "..." << cgiter << " ( " << omega1 << " )";
    }

    if (omega1 < epsilon2 * omega_z) {
      optimality = 1;
      break;
    }

    scale = omega1 * inv_omega2;
    dscal_(&n, &scale, p, &inc);
    daxpy_(&n, &one, r, &inc, p, &inc);
    omega_p = ddot_(&n, p, &inc, p, &inc);
  }
  if (VERBOSE_CGLS) {
    cout << "...Done." << endl;
  }
  tictoc.stop();
  if (VERB > 4) {
    cerr << "CGLS converged in " << cgiter << " iteration(s) and "
        << tictoc.getCPUTimeStr() << " CPU seconds." << endl;
  }
  return optimality;
}

int L2_SVM_MFN(const AlgIn& data, const options& Options,
               std::vector<double>& Weights,
               std::vector<double>& Outputs, double cpos, double cneg) {
  /* Disassemble the structures */
  Timer tictoc;
  const double* Y = data.Y.data();
  int n = static_cast<int>(Weights.size());
  const int m = data.m;
  double lambda = Options.lambda;
  double epsilon = BIG_EPSILON;
  int cgitermax = SMALL_CGITERMAX;
  double* w = Weights.data();
  double* o = Outputs.data();
  double F_old = 0.0;
  double F = 0.0;
  double diff = 0.0;
  int ini = 0;
  int n0 = n-1;
  int inc = 1;
  // active examples first, inactive ones from the back
  std::vector<int> activeSubset(static_cast<size_t>(m));
  int* activeSet = activeSubset.data();
  // initialize
  F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
  int active = 0;
  int inactive = m - 1; // l-1
  for (int i = 0; i < m; i++) {
    diff = 1 - Y[i] * o[i];
    if (diff > 0) {
      activeSet[active] = i;
      active++;
      // C[i]
      F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
    } else {
      activeSet[inactive] = i;
      inactive--;
    }
  }
  int iter = 0;
  int opt = 0;
  int opt2 = 0;
  std::vector<double> weightsBar(static_cast<size_t>(n));
  std::vector<double> outputsBar(static_cast<size_t>(m));
  double* w_bar = weightsBar.data();
  double* o_bar = outputsBar.data();
  double delta = 0.0;
  int ii = 0;
  SolverWorkspace workspace(m, n);
  while (iter < Options.mfnitermax) {
    iter++;
    if (VERB > 4) {
      cerr << "L2_SVM_MFN Iteration# " << iter << " (" << active
          << " active examples, " << " objective_value = " << F << ")"
          << endl;
    }
    memcpy(w_bar, w, sizeof(double)*static_cast<std::size_t>(n));
    memcpy(o_bar, o, sizeof(double)*static_cast<std::size_t>(m));
    opt = Options.directSolve
              ? directLeastSquares(data, lambda, activeSet, active, w_bar,
                                   o_bar, cpos, cneg, workspace)
              : CGLS(data, lambda, cgitermax, epsilon, activeSet, active,
                     w_bar, o_bar, cpos, cneg, workspace);
    // outputs of the inactive examples, which CGLS did not update
    scoreRows(activeSet + active, m - active, data, n0, w_bar, o_bar);
    if (ini == 0) {
      cgitermax = CGITERMAX;
      ini = 1;
    };
    opt2 = 1;
    for (int i = 0; i < m; i++) {
      ii = activeSet[i];
      if (i < active) {
        opt2 = (opt2 && (Y[ii] * o_bar[ii] <= 1 + epsilon));
      } else {
        opt2 = (opt2 && (Y[ii] * o_bar[ii] >= 1 - epsilon));
      }
      if (opt2 == 0) {
        break;
      }
    }
    if (opt && opt2) { // l
      if (epsilon == BIG_EPSILON) {
        epsilon = Options.epsilon;
        if (VERB > 4) {
          cerr << "  epsilon = " << BIG_EPSILON
            << " case converged (speedup heuristic 2). Continuing with epsilon="
            << EPSILON << endl;
        }
        continue;
      } else {
        memcpy(w, w_bar, sizeof(double)*static_cast<std::size_t>(n));
        memcpy(o, o_bar, sizeof(double)*static_cast<std::size_t>(m));
        if (VERB > 3) {
          tictoc.stop();
          cerr << "L2_SVM_MFN converged (optimality) in " << iter
              << " iteration(s) and " << tictoc.getCPUTimeStr() << " CPU seconds. \n"
              << endl;
        }
        return 1;
      }
    }
    delta = line_search(w, w_bar, lambda, o, o_bar, Y, n, m, cpos, cneg,
                        workspace.deltas.data());
    F_old = F;
    double delta2 = 1-delta;
    dscal_(&n, &delta2, w, &inc);
    daxpy_(&n, &delta, w_bar, &inc, w, &inc);
    F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
    active = 0;
    inactive = m - 1;
    for (int i = 0; i < m; i++) {
      o[i] += delta * (o_bar[i] - o[i]);
      diff = 1 - Y[i] * o[i];
      if (diff > 0) {
        activeSet[active] = i;
        active++;
        F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
      } else {
        activeSet[inactive] = i;
        inactive--;
      }
    }
      if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
      return 2;
    }
  }
  return 0;
}

double line_search(double* w, double* w_bar, double lambda, double* o,
                   double* o_bar, const double* Y, int d, /* data dimensionality -- 'n' */
                   int l, double cpos, double cneg, Delta* deltas){
  int i = 0;
  double omegaL = 0.0;
  double omegaR = 0.0;
  double diff = 0.0;
  for (int i = d; i--;) {
    diff = w_bar[i] - w[i];
    omegaL += w[i] * diff;
    omegaR += w_bar[i] * diff;
  }
  omegaL = lambda * omegaL;
  omegaR = lambda * omegaR;
  double L = omegaL;
  double R = omegaR;
  int ii = 0;
  double d2 = 0.0;

  int p = 0;
  for (i = 0; i < l; i++) {
    diff = Y[i] * (o_bar[i] - o[i]);
    if (Y[i] * o[i] < 1) {
      d2 = ((Y[i]==1)? cpos : cneg) * (o_bar[i] - o[i]);
      L += (o[i] - Y[i]) * d2;
      R += (o_bar[i] - Y[i]) * d2;
      if (diff > 0) {
        deltas[p].delta = (1 - Y[i] * o[i]) / diff;
        deltas[p].index = i;
        deltas[p].s = -1;
        p++;
      }
    } else {
      if (diff < 0) {
        deltas[p].delta = (1 - Y[i] * o[i]) / diff;
        deltas[p].index = i;
        deltas[p].s = 1;
        p++;
      }
    }
  }
  sort(deltas, deltas + p);
  double delta_prime = 0.0;
  for (i = 0; i < p; i++) {
    delta_prime = L + deltas[i].delta * (R - L);
    if (delta_prime >= 0) {
      break;
    }
    ii = deltas[i].index;
    diff = (deltas[i].s) * ((Y[ii]==1)? cpos : cneg) * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  return (-L / (R - L));
}

