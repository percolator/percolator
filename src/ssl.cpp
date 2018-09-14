/*******************************************************************************
 * SVMlin
 * Copyright (c) 2006 Vikas Sindhwani at the University of Chicago.
 * Adapted to Percolator by Lukas Käll at the University of Washington
 *
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

extern "C" {
  extern double dnrm2_(int *, double *, int *);
  extern double ddot_(int *, double *, int *, double *, int *);
  extern int daxpy_(int *, double *, double *, int *, double *, int *);
  extern int dscal_(int *, double *, double *, int *);
  extern int dgemv_(char *, int *, int *,
		    double *, double *, int *,
		    double *, int *,  double *,
		    double *, int *);
}

#define VERBOSE 1
#define LOG2(x) 1.4426950408889634*log(x)
// for compatibility issues, not using log2

AlgIn::AlgIn(const int size, const int numFeat) {
  vals = new double*[size];
  Y = new double[size];
  C = new double[size];
  n = numFeat;
  positives = 0;
  negatives = 0;
}
AlgIn::~AlgIn() {
  delete[] vals;
  delete[] Y;
  delete[] C;
}

double cglsFun1(int active, int* J, const double* Y,
		double* set2, int n, double* q, 
		double* p, double cpos, double cneg){
  double omega_q = 0.0;
  int inc = 1;
  int i = 0;

  char trans = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  dgemv_(&trans, &n, &active,
	 &alpha, set2, &n,
	 p, &inc, &beta, q, &inc);

  for (i = 0; i < active; i++) {
    omega_q += ((Y[J[i]]==1)? cpos : cneg) * (q[i]) * (q[i]);
  }

  return(omega_q);
}

void cglsFun2(int active, int* J, const double* Y,
	      double* set2, int n0, int n, double* q, 
	      double* o, double* z, double* r, 
	      double cpos, double cneg){
  int i;
  int inc = 1;
  int ind = 0;
  double temp = 0.0;
  
  for (i = 0; i < active; i++) {
    o[J[i]] += q[i];
    z[i] -= ((Y[J[i]]==1)? cpos : cneg) * q[i];
    daxpy_(&n, &(z[i]), set2 + i * n, &inc, r, &inc);
  }
}

int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
         const double epsilon, const struct vector_int* Subset,
         struct vector_double* Weights, struct vector_double* Outputs,
         double cpos, double cneg) {
  if (VERBOSE_CGLS) {
    cout << "CGLS starting..." << endl;
  }
  /* Disassemble the structures */
  timer tictoc;
  tictoc.restart();
  int active = Subset->d;
  int* J = Subset->vec;
  double** set = data.vals;
  const double* Y = data.Y;
  //const double* C = data.C;
  int n = data.n;
  //  int m  = pSet->size();
  double* beta = Weights->vec;
  double* o = Outputs->vec;
  // initialize z
  double* z = new double[active];
  double* q = new double[active];
  int ii = 0;
  register int i, j;
  int n0 = n-1;
  int inc = 1;
  double one = 1;
  double negLambda = -lambda;
  int rowStart = 0;
  double* set2 = new double[n*active];
  double* r = new double[n];
  for (i = n; i--;) {
    r[i] = 0.0;
  }
  for (i = 0; i < active; i++) {
    ii = J[i];
    z[i] = ((Y[ii]==1)? cpos : cneg) * (Y[ii] - o[ii]);
    rowStart = i * n;
    memcpy(set2 + rowStart, set[ii], sizeof(double)*n0);
    set2[rowStart + n0] = 1.0;
    daxpy_(&n, &(z[i]), set2 + i*n, &inc, r, &inc);
  }
  double* p = new double[n];
  daxpy_(&n, &negLambda, beta, &inc, r, &inc);
  memcpy(p, r, sizeof(double)*n);
  double omega1 = ddot_(&n, r, &inc, r, &inc);
  // for (i = active; i--;) {
  //   ii = J[i];
  //   // C[ii]
  //   z[i] = ((Y[ii]==1)? cpos : cneg) * (Y[ii] - o[ii]);
  // }
  // double* r = new double[n];
  // for (i = n; i--;) {
  //   r[i] = 0.0;
  // }
  // for (j = 0; j < active; j++) {
  //   const double* val = set[J[j]];
  //   for (i = n - 1; i--;) {
  //     r[i] += val[i] * z[j];
  //   }
  //   r[n - 1] += z[j];
  // }
  // double* p = new double[n];
  // double omega1 = 0.0;
  // for (i = n; i--;) {
  //   r[i] -= lambda * beta[i];
  //   p[i] = r[i];
  //   omega1 += r[i] * r[i];
  // }
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
    omega_q = cglsFun1(active, J, Y, set2, n, q, p, cpos, cneg);
    gamma = omega1 / (lambda * omega_p + omega_q);
    inv_omega2 = 1 / omega1;

    memcpy(r, beta, sizeof(double)*n);
    dscal_(&n, &negLambda, r, &inc);

    daxpy_(&n, &gamma, p, &inc, beta, &inc);
    dscal_(&active, &gamma, q, &inc);

    cglsFun2(active, J, Y, set2,
	     n0, n, q, o, z, r, cpos, cneg);

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

    // omega_q = 0.0;
    // double t = 0.0;
    // //    register int i,j;
    // // #pragma omp parallel for private(i,j)
    // for (i = 0; i < active; i++) {
    //   ii = J[i];
    //   t = 0.0;
    //   const double* val = set[ii];
    //   for (j = 0; j < n - 1; j++) {
    //     t += val[j] * p[j];
    //   }
    //   t += p[n - 1];
    //   q[i] = t;
    //   // C[ii]
    //   omega_q += ((Y[ii]==1)? cpos : cneg) * t * t;
    // }
    // gamma = omega1 / (lambda * omega_p + omega_q);
    // inv_omega2 = 1 / omega1;
    // for (int i = n; i--;) {
    //   r[i] = 0.0;
    //   beta[i] += gamma * p[i];
    // }
    // omega_z = 0.0;
    // for (int i = active; i--;) {
    //   ii = J[i];
    //   o[ii] += gamma * q[i];
    //   // C[ii]
    //   z[i] -= gamma * ((Y[ii]==1)? cpos : cneg) * q[i];
    //   omega_z += z[i] * z[i];
    // }
    // for (register int j = 0; j < active; j++) {
    //   ii = J[j];
    //   t = z[j];
    //   const double* val = set[ii];
    //   for (register int i = 0; i < n - 1; i++) {
    //     r[i] += val[i] * t;
    //   }
    //   r[n - 1] += t;
    // }
    // omega1 = 0.0;
    // for (int i = n; i--;) {
    //   r[i] -= lambda * beta[i];
    //   omega1 += r[i] * r[i];
    // }
    // if (VERBOSE_CGLS) {
    //   cout << "..." << cgiter << " ( " << omega1 << " )";
    // }
    // if (omega1 < epsilon2 * omega_z) {
    //   optimality = 1;
    //   break;
    // }
    // omega_p = 0.0;
    // scale = omega1 * inv_omega2;
    // for (int i = n; i--;) {
    //   p[i] = r[i] + p[i] * scale;
    //   omega_p += p[i] * p[i];
    // }
  }
  if (VERBOSE_CGLS) {
    cout << "...Done." << endl;
  }
  tictoc.stop();
  if (VERB > 4) {
    cerr << "CGLS converged in " << cgiter << " iteration(s) and "
        << tictoc.time() << " seconds." << endl;
  }
  delete[] z;
  delete[] q;
  delete[] r;
  delete[] p;
  delete[] set2;
  return optimality;
}

int L2_SVM_MFN(const AlgIn& data, struct options* Options,
               struct vector_double* Weights,
               struct vector_double* Outputs, double cpos, double cneg) {
  /* Disassemble the structures */
  timer tictoc;
  tictoc.restart();
  double** set = data.vals;
  const double* Y = data.Y;
  //const double* C = data.C;
  int n = Weights->d;
  const int m = data.m;
  double lambda = Options->lambda;
  double epsilon = BIG_EPSILON;
  int cgitermax = SMALL_CGITERMAX;
  double* w = Weights->vec;
  double* o = Outputs->vec;
  double F_old = 0.0;
  double F = 0.0;
  double diff = 0.0;
  int ini = 0;
  int n0 = n-1;
  int inc = 1;
  vector_int* ActiveSubset = new vector_int[1];
  ActiveSubset->vec = new int[m];
  ActiveSubset->d = m;
  // initialize
  F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
  // for (int i = 0; i < n; i++) {
  //   F += w[i] * w[i];
  // }
  // F = 0.5 * lambda * F;
  int active = 0;
  int inactive = m - 1; // l-1
  for (int i = 0; i < m; i++) {
    diff = 1 - Y[i] * o[i];
    if (diff > 0) {
      ActiveSubset->vec[active] = i;
      active++;
      // C[i]
      F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
    } else {
      ActiveSubset->vec[inactive] = i;
      inactive--;
    }
  }
  ActiveSubset->d = active;
  int iter = 0;
  int opt = 0;
  int opt2 = 0;
  vector_double* Weights_bar = new vector_double[1];
  vector_double* Outputs_bar = new vector_double[1];
  double* w_bar = new double[n];
  double* o_bar = new double[m];
  Weights_bar->vec = w_bar;
  Outputs_bar->vec = o_bar;
  Weights_bar->d = n;
  Outputs_bar->d = m;
  double delta = 0.0;
  double t = 0.0;
  int ii = 0;
  while (iter < Options->mfnitermax) {
    iter++;
    if (VERB > 4) {
      cerr << "L2_SVM_MFN Iteration# " << iter << " (" << active
          << " active examples, " << " objective_value = " << F << ")"
          << endl;
    }
    memcpy(w_bar, w, sizeof(double)*n);
    memcpy(o_bar, o, sizeof(double)*m);
    // for (int i = n; i--;) {
    //   w_bar[i] = w[i];
    // }
    // for (int i = m; i--;) {
    //   o_bar[i] = o[i];
    // }
    opt = CGLS(data,
               lambda,
               cgitermax,
               epsilon,
               ActiveSubset,
               Weights_bar,
               Outputs_bar, cpos, cneg);
    for (register int i = active; i < m; i++) {
      ii = ActiveSubset->vec[i];
      o_bar[ii] = ddot_(&n0, set[ii], &inc, w_bar, &inc) + w_bar[n - 1];
      // const double* val = set[ii];
      // t = w_bar[n - 1];
      // for (register int j = n - 1; j--;) {
      //   t += val[j] * w_bar[j];
      // }
      // o_bar[ii] = t;
    }
    if (ini == 0) {
      cgitermax = CGITERMAX;
      ini = 1;
    };
    opt2 = 1;
    for (int i = 0; i < m; i++) {
      ii = ActiveSubset->vec[i];
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
        epsilon = Options->epsilon;
        if (VERB > 4) {
          cerr << "  epsilon = " << BIG_EPSILON
            << " case converged (speedup heuristic 2). Continuing with epsilon="
            << EPSILON << endl;
        }
        continue;
      } else {
	memcpy(w, w_bar, sizeof(double)*n);
	memcpy(o, o_bar, sizeof(double)*m);
        // for (int i = n; i--;) {
        //   w[i] = w_bar[i];
        // }
        // for (int i = m; i--;) {
        //   o[i] = o_bar[i];
        // }
        delete[] ActiveSubset->vec;
        delete[] ActiveSubset;
        delete[] o_bar;
        delete[] w_bar;
        delete[] Weights_bar;
        delete[] Outputs_bar;
        tictoc.stop();
        if (VERB > 3) {
          cerr << "L2_SVM_MFN converged (optimality) in " << iter
              << " iteration(s) and " << tictoc.time() << " seconds. \n"
              << endl;
        }
        return 1;
      }
    }
    delta = line_search_nrOne(w, w_bar, lambda, o, o_bar, Y, n, m, cpos, cneg);
    // delta = line_search(w, w_bar, lambda, o, o_bar, Y, NULL, n, m, cpos, cneg);
    F_old = F;
    double delta2 = 1-delta;
    dscal_(&n, &delta2, w, &inc);
    daxpy_(&n, &delta, w_bar, &inc, w, &inc);
    F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
    // F = 0.0;
    // for (int i = n; i--;) {
    //   w[i] += delta * (w_bar[i] - w[i]);
    //   F += w[i] * w[i];
    // }
    // F = 0.5 * lambda * F;
    active = 0;
    inactive = m - 1;
    for (int i = 0; i < m; i++) {
      o[i] += delta * (o_bar[i] - o[i]);
      diff = 1 - Y[i] * o[i];
      if (diff > 0) {
        ActiveSubset->vec[active] = i;
        active++;
        F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
      } else {
        ActiveSubset->vec[inactive] = i;
        inactive--;
      }
    }
    ActiveSubset->d = active;
    if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
      //    cout << "L2_SVM_MFN converged (rel. criterion) in " << iter << " iterations and "<< tictoc.time() << " seconds. \n" << endl;
      return 2;
    }
  }
  delete[] ActiveSubset->vec;
  delete[] ActiveSubset;
  delete[] o_bar;
  delete[] w_bar;
  delete[] Weights_bar;
  delete[] Outputs_bar;
  tictoc.stop();
  //  cout << "L2_SVM_MFN converged (max iter exceeded) in " << iter << " iterations and "<< tictoc.time() << " seconds. \n" << endl;
  return 0;
}

double line_search_nrOne(double* w, double* w_bar, double lambda, double* o,
			 double* o_bar, const double* Y, int d, /* data dimensionality -- 'n' */
			 int l, double cpos, double cneg){
		   // ) { /* number of examples */
  int inc = 1;
  int i = 0;
  double omegaL = 0.0; // lambda * ddot_(&d, w, &inc, w_barD, &inc);
  double omegaR = 0.0; // lambda * ddot_(&d, w_bar, &inc, w_barD, &inc);
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

  Delta* deltas = new Delta[l];
  int p = 0;
  for (i = 0; i < l; i++) {
    diff = Y[i] * (o_bar[i] - o[i]);
    if (Y[i] * o[i] < 1) {
      // d2 = C[i] * (o_bar[i] - o[i]);
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
    // diff = (deltas[i].s) * C[ii] * (o_bar[ii] - o[ii]);
    diff = (deltas[i].s) * ((Y[ii]==1)? cpos : cneg) * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  delete[] deltas;
  return (-L / (R - L));
}

double line_search(double* w, double* w_bar, double lambda, double* o,
                   double* o_bar, const double* Y, const double* C_, int d, /* data dimensionality -- 'n' */
                   int l, double cpos, double cneg) { /* number of examples */
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
  double L = 0.0;
  double R = 0.0;
  int ii = 0;
  for (int i = 0; i < l; i++) {
    if (Y[i] * o[i] < 1) {
      // C[i]
      diff = ((Y[i]==1)? cpos : cneg) * (o_bar[i] - o[i]);
      L += (o[i] - Y[i]) * diff;
      R += (o_bar[i] - Y[i]) * diff;
    }
  }
  L += omegaL;
  R += omegaR;
  Delta* deltas = new Delta[l];
  int p = 0;
  for (int i = 0; i < l; i++) {
    diff = Y[i] * (o_bar[i] - o[i]);
    if (Y[i] * o[i] < 1) {
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
  for (int i = 0; i < p; i++) {
    delta_prime = L + deltas[i].delta * (R - L);
    if (delta_prime >= 0) {
      break;
    }
    ii = deltas[i].index;
    // C[ii]
    diff = (deltas[i].s) * ((Y[ii]==1)? cpos : cneg) * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  delete[] deltas;
  return (-L / (R - L));
}
/********************** UTILITIES ********************/

void Clear(struct data* a) {
  delete[] a->val;
  delete[] a->rowptr;
  delete[] a->colind;
  delete[] a->Y;
  delete[] a->C;
  delete a;
  return;
}
void Clear(struct vector_double* c) {
  delete[] c->vec;
  delete[] c;
  return;
}
void Clear(struct vector_int* c) {
  delete[] c->vec;
  delete[] c;
  return;
}
void Clear(struct options* opt) {
  delete[] opt;
  delete[] opt;
  return;
}
