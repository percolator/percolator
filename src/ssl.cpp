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
#include <cstring>
#include <set>
#include <vector>
#include <ctype.h>
using namespace std;
#include "Globals.h"
#include "ssl.h"

#include <stdarg.h>
#include <omp.h>

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

#define THREADS (Globals::getInstance()->getNumThreads())
#define FUN2
#define LFUN2
// for compatibility issues, not using log2

/********************** Parallel reduction class, based off liblinear class
Author: John T. Halloran
Affiliation: UC Davis
Date: May 2017
********************/
class Reduce_Vectors
{
 public:
  Reduce_Vectors(int size);
  ~Reduce_Vectors();

  void init(void);
  void sum_scale_x(double scalar, double* x);
  void reduce_sum(double* v);
  void reduce_sum_cgls(double* v);

 private:
  int nr_thread;
  int size;
  int size0;
  double **tmp_array;
};

Reduce_Vectors::Reduce_Vectors(int size)
{
  #ifdef NOMAP
  nr_thread = omp_get_max_threads();
  #else
  nr_thread = 3 * THREADS;
  #endif
  this->size = size;
  this->size0 = size-1;
  tmp_array = new double*[nr_thread];
  for(int i = 0; i < nr_thread; i++)
    tmp_array[i] = new double[size];
}

Reduce_Vectors::~Reduce_Vectors(void)
{
  for(int i = 0; i < nr_thread; i++)
    delete[] tmp_array[i];
  delete[] tmp_array;
}

void Reduce_Vectors::init(void)
{
#pragma omp parallel for schedule(static)
  for(int i = 0; i < size; i++)
    for(int j = 0; j < nr_thread; j++)
      tmp_array[j][i] = 0.0;
}

void Reduce_Vectors::sum_scale_x(double scalar, double* x)
{
  int thread_id = omp_get_thread_num();
  int inc = 1;

// #ifdef BLASC
//   cblas_daxpy(size0, scalar, x, inc, tmp_array[thread_id], inc);
// #else
  daxpy_(&size, &scalar, x, &inc, tmp_array[thread_id], &inc);
  // #endif
  // tmp_array[thread_id][size-1] += scalar;

// #ifdef BLASC
//   cblas_daxpy(size, scalar, x, inc, tmp_array[thread_id], inc);
// #else
//   daxpy_(&size, &scalar, x, &inc, tmp_array[thread_id], &inc);
// #endif
}

void Reduce_Vectors::reduce_sum(double* v)
{
#pragma omp parallel for schedule(static)
  for(int i = 0; i < size; i++)
    {
      v[i] = 0;
      for(int j = 0; j < nr_thread; j++)
	v[i] += tmp_array[j][i];
    }
}

void Reduce_Vectors::reduce_sum_cgls(double* v)
{
#pragma omp parallel for schedule(static)
  for(int i = 0; i < size; i++)
    {
      for(int j = 0; j < nr_thread; j++)
	v[i] += tmp_array[j][i];
    }
}
/////////////////////////

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


// //////////////////////////////////////////////////////////////////////////
// ////////////////  John: leaving original, unmodified L2-SVM-MFN code here
// ////////////////        for completeness
// //////////////////////////////////////////////////////////////////////////
// int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
//          const double epsilon, const struct vector_int* Subset,
//          struct vector_double* Weights, struct vector_double* Outputs) {
//   if (VERBOSE_CGLS) {
//     cout << "CGLS starting..." << endl;
//   }
//   /* Disassemble the structures */
//   timer tictoc;
//   tictoc.restart();
//   int active = Subset->d;
//   int* J = Subset->vec;
//   const double** set = data.vals;
//   const double* Y = data.Y;
//   const double* C = data.C;
//   const int n = data.n;
//   //  int m  = pSet->size();
//   double* beta = Weights->vec;
//   double* o = Outputs->vec;
//   // initialize z
//   double* z = new double[active];
//   double* q = new double[active];
//   int ii = 0;
//   register int i, j;
//   for (i = active; i--;) {
//     ii = J[i];
//     z[i] = C[ii] * (Y[ii] - o[ii]);
//   }
//   double* r = new double[n];
//   for (i = n; i--;) {
//     r[i] = 0.0;
//   }
//   for (j = 0; j < active; j++) {
//     const double* val = set[J[j]];
//     for (i = n - 1; i--;) {
//       r[i] += val[i] * z[j];
//     }
//     r[n - 1] += z[j];
//   }
//   double* p = new double[n];
//   double omega1 = 0.0;
//   for (i = n; i--;) {
//     r[i] -= lambda * beta[i];
//     p[i] = r[i];
//     omega1 += r[i] * r[i];
//   }
//   double omega_p = omega1;
//   double omega_q = 0.0;
//   double inv_omega2 = 1 / omega1;
//   double scale = 0.0;
//   double omega_z = 0.0;
//   double gamma = 0.0;
//   int cgiter = 0;
//   int optimality = 0;
//   double epsilon2 = epsilon * epsilon;
//   // iterate
//   while (cgiter < cgitermax) {
//     cgiter++;
//     omega_q = 0.0;
//     double t = 0.0;
//     //    register int i,j;
//     // #pragma omp parallel for private(i,j)
//     for (i = 0; i < active; i++) {
//       ii = J[i];
//       t = 0.0;
//       const double* val = set[ii];
//       for (j = 0; j < n - 1; j++) {
//         t += val[j] * p[j];
//       }
//       t += p[n - 1];
//       q[i] = t;
//       omega_q += C[ii] * t * t;
//     }
//     gamma = omega1 / (lambda * omega_p + omega_q);
//     inv_omega2 = 1 / omega1;
//     for (int i = n; i--;) {
//       r[i] = 0.0;
//       beta[i] += gamma * p[i];
//     }
//     omega_z = 0.0;
//     for (int i = active; i--;) {
//       ii = J[i];
//       o[ii] += gamma * q[i];
//       z[i] -= gamma * C[ii] * q[i];
//       omega_z += z[i] * z[i];
//     }
//     for (register int j = 0; j < active; j++) {
//       ii = J[j];
//       t = z[j];
//       const double* val = set[ii];
//       for (register int i = 0; i < n - 1; i++) {
//         r[i] += val[i] * t;
//       }
//       r[n - 1] += t;
//     }
//     omega1 = 0.0;
//     for (int i = n; i--;) {
//       r[i] -= lambda * beta[i];
//       omega1 += r[i] * r[i];
//     }
//     if (VERBOSE_CGLS) {
//       cout << "..." << cgiter << " ( " << omega1 << " )";
//     }
//     if (omega1 < epsilon2 * omega_z) {
//       optimality = 1;
//       break;
//     }
//     omega_p = 0.0;
//     scale = omega1 * inv_omega2;
//     for (int i = n; i--;) {
//       p[i] = r[i] + p[i] * scale;
//       omega_p += p[i] * p[i];
//     }
//   }
//   if (VERBOSE_CGLS) {
//     cout << "...Done." << endl;
//   }
//   tictoc.stop();
//   if (VERB > 4) {
//     cerr << "CGLS converged in " << cgiter << " iteration(s) and "
//         << tictoc.time() << " seconds." << endl;
//   }
//   delete[] z;
//   delete[] q;
//   delete[] r;
//   delete[] p;
//   return optimality;
// }

// int L2_SVM_MFN(const AlgIn& data, struct options* Options,
//                struct vector_double* Weights,
//                struct vector_double* Outputs) {
//   /* Disassemble the structures */
//   timer tictoc;
//   tictoc.restart();
//   const double** set = data.vals;
//   const double* Y = data.Y;
//   const double* C = data.C;
//   const int n = Weights->d;
//   const int m = data.m;
//   double lambda = Options->lambda;
//   double epsilon = BIG_EPSILON;
//   int cgitermax = SMALL_CGITERMAX;
//   double* w = Weights->vec;
//   double* o = Outputs->vec;
//   double F_old = 0.0;
//   double F = 0.0;
//   double diff = 0.0;
//   int ini = 0;
//   vector_int* ActiveSubset = new vector_int[1];
//   ActiveSubset->vec = new int[m];
//   ActiveSubset->d = m;
//   // initialize
//   for (int i = 0; i < n; i++) {
//     F += w[i] * w[i];
//   }
//   F = 0.5 * lambda * F;
//   int active = 0;
//   int inactive = m - 1; // l-1
//   for (int i = 0; i < m; i++) {
//     diff = 1 - Y[i] * o[i];
//     if (diff > 0) {
//       ActiveSubset->vec[active] = i;
//       active++;
//       F += 0.5 * C[i] * diff * diff;
//     } else {
//       ActiveSubset->vec[inactive] = i;
//       inactive--;
//     }
//   }
//   ActiveSubset->d = active;
//   int iter = 0;
//   int opt = 0;
//   int opt2 = 0;
//   vector_double* Weights_bar = new vector_double[1];
//   vector_double* Outputs_bar = new vector_double[1];
//   double* w_bar = new double[n];
//   double* o_bar = new double[m];
//   Weights_bar->vec = w_bar;
//   Outputs_bar->vec = o_bar;
//   Weights_bar->d = n;
//   Outputs_bar->d = m;
//   double delta = 0.0;
//   double t = 0.0;
//   int ii = 0;
//   while (iter < Options->mfnitermax) {
//     iter++;
//     if (VERB > 4) {
//       cerr << "L2_SVM_MFN Iteration# " << iter << " (" << active
//           << " active examples, " << " objective_value = " << F << ")"
//           << endl;
//     }
//     for (int i = n; i--;) {
//       w_bar[i] = w[i];
//     }
//     for (int i = m; i--;) {
//       o_bar[i] = o[i];
//     }
//     opt = CGLS(data,
//                lambda,
//                cgitermax,
//                epsilon,
//                ActiveSubset,
//                Weights_bar,
//                Outputs_bar);
//     for (register int i = active; i < m; i++) {
//       ii = ActiveSubset->vec[i];
//       const double* val = set[ii];
//       t = w_bar[n - 1];
//       for (register int j = n - 1; j--;) {
//         t += val[j] * w_bar[j];
//       }
//       o_bar[ii] = t;
//     }
//     if (ini == 0) {
//       cgitermax = CGITERMAX;
//       ini = 1;
//     };
//     opt2 = 1;
//     for (int i = 0; i < m; i++) {
//       ii = ActiveSubset->vec[i];
//       if (i < active) {
//         opt2 = (opt2 && (Y[ii] * o_bar[ii] <= 1 + epsilon));
//       } else {
//         opt2 = (opt2 && (Y[ii] * o_bar[ii] >= 1 - epsilon));
//       }
//       if (opt2 == 0) {
//         break;
//       }
//     }
//     if (opt && opt2) { // l
//       if (epsilon == BIG_EPSILON) {
//         epsilon = Options->epsilon;
//         if (VERB > 4) {
//           cerr << "  epsilon = " << BIG_EPSILON
//             << " case converged (speedup heuristic 2). Continuing with epsilon="
//             << EPSILON << endl;
//         }
//         continue;
//       } else {
//         for (int i = n; i--;) {
//           w[i] = w_bar[i];
//         }
//         for (int i = m; i--;) {
//           o[i] = o_bar[i];
//         }
//         delete[] ActiveSubset->vec;
//         delete[] ActiveSubset;
//         delete[] o_bar;
//         delete[] w_bar;
//         delete[] Weights_bar;
//         delete[] Outputs_bar;
//         tictoc.stop();
//         if (VERB > 3) {
//           cerr << "L2_SVM_MFN converged (optimality) in " << iter
//               << " iteration(s) and " << tictoc.time() << " seconds. \n"
//               << endl;
//         }
//         return 1;
//       }
//     }
//     delta = line_search(w, w_bar, lambda, o, o_bar, Y, C, n, m);
//     F_old = F;
//     F = 0.0;
//     for (int i = n; i--;) {
//       w[i] += delta * (w_bar[i] - w[i]);
//       F += w[i] * w[i];
//     }
//     F = 0.5 * lambda * F;
//     active = 0;
//     inactive = m - 1;
//     for (int i = 0; i < m; i++) {
//       o[i] += delta * (o_bar[i] - o[i]);
//       diff = 1 - Y[i] * o[i];
//       if (diff > 0) {
//         ActiveSubset->vec[active] = i;
//         active++;
//         F += 0.5 * C[i] * diff * diff;
//       } else {
//         ActiveSubset->vec[inactive] = i;
//         inactive--;
//       }
//     }
//     ActiveSubset->d = active;
//     if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
//       //    cout << "L2_SVM_MFN converged (rel. criterion) in " << iter << " iterations and "<< tictoc.time() << " seconds. \n" << endl;
//       return 2;
//     }
//   }
//   delete[] ActiveSubset->vec;
//   delete[] ActiveSubset;
//   delete[] o_bar;
//   delete[] w_bar;
//   delete[] Weights_bar;
//   delete[] Outputs_bar;
//   tictoc.stop();
//   //  cout << "L2_SVM_MFN converged (max iter exceeded) in " << iter << " iterations and "<< tictoc.time() << " seconds. \n" << endl;
//   return 0;
// }

// double line_search(double* w, double* w_bar, double lambda, double* o,
//                    double* o_bar, const double* Y, const double* C, int d, /* data dimensionality -- 'n' */
//                    int l) { /* number of examples */
//   double omegaL = 0.0;
//   double omegaR = 0.0;
//   double diff = 0.0;
//   for (int i = d; i--;) {
//     diff = w_bar[i] - w[i];
//     omegaL += w[i] * diff;
//     omegaR += w_bar[i] * diff;
//   }
//   omegaL = lambda * omegaL;
//   omegaR = lambda * omegaR;
//   double L = 0.0;
//   double R = 0.0;
//   int ii = 0;
//   for (int i = 0; i < l; i++) {
//     if (Y[i] * o[i] < 1) {
//       diff = C[i] * (o_bar[i] - o[i]);
//       L += (o[i] - Y[i]) * diff;
//       R += (o_bar[i] - Y[i]) * diff;
//     }
//   }
//   L += omegaL;
//   R += omegaR;
//   Delta* deltas = new Delta[l];
//   int p = 0;
//   for (int i = 0; i < l; i++) {
//     diff = Y[i] * (o_bar[i] - o[i]);
//     if (Y[i] * o[i] < 1) {
//       if (diff > 0) {
//         deltas[p].delta = (1 - Y[i] * o[i]) / diff;
//         deltas[p].index = i;
//         deltas[p].s = -1;
//         p++;
//       }
//     } else {
//       if (diff < 0) {
//         deltas[p].delta = (1 - Y[i] * o[i]) / diff;
//         deltas[p].index = i;
//         deltas[p].s = 1;
//         p++;
//       }
//     }
//   }
//   sort(deltas, deltas + p);
//   double delta_prime = 0.0;
//   for (int i = 0; i < p; i++) {
//     delta_prime = L + deltas[i].delta * (R - L);
//     if (delta_prime >= 0) {
//       break;
//     }
//     ii = deltas[i].index;
//     diff = (deltas[i].s) * C[ii] * (o_bar[ii] - o[ii]);
//     L += diff * (o[ii] - Y[ii]);
//     R += diff * (o_bar[ii] - Y[ii]);
//   }
//   delete[] deltas;
//   return (-L / (R - L));
// }

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

///////////////////  Optimized SVM solvers, for both single and multi-threading

/********************** L2-SVM-MFN optimized implementations
Author: John T. Halloran
Affiliation: UC Davis
Date: May 2017
********************/
double cglsFun1_nrOne(int active, int* J, const double* C, 
		      double* set2, int n, double* q, 
		      double* p){
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
    omega_q += C[J[i]] * (q[i]) * (q[i]);
  }

  return(omega_q);
}

void cglsFun2_nrOne(int active, int* J, const double* C, 
			   double* set2, int n0, int n, double* q, 
			   double* o, double* z, double* r){
  int i;
  int inc = 1;
  int ind = 0;
  double temp = 0.0;
  
  for (i = 0; i < active; i++) {
    o[J[i]] += q[i];
    z[i] -= C[J[i]] * q[i];
    daxpy_(&n, &(z[i]), set2 + i * n, &inc, r, &inc);
  }
}

int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
	 const double epsilon, const struct vector_int* Subset,
	 struct vector_double* Weights, struct vector_double* Outputs) {
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
  const double* C = data.C;
  // double* tempZ = new double[active];
  int n = data.n;
  //  int m  = pSet->size();
  double* beta = Weights->vec;
  double* o = Outputs->vec;
  // initialize z
  double* z = new double[active];
  double* q = new double[active];
  int ii = 0;
  int i, j;
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
    z[i] = C[ii] * (Y[ii] - o[ii]);
    rowStart = i * n;
    memcpy(set2 + rowStart, set[ii], sizeof(double)*n0);
    set2[rowStart + n0] = 1.0;
    daxpy_(&n, &(z[i]), set2 + i*n, &inc, r, &inc);
  }
  // for (j = 0; j < active; j++) {
  //   daxpy_(&n, &(z[j]), set2 + j*n, &inc, r, &inc);
  // }
  double* p = new double[n];
  daxpy_(&n, &negLambda, beta, &inc, r, &inc);
  memcpy(p, r, sizeof(double)*n);
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
    omega_q = cglsFun1_nrOne(active, J, C, set2, n, q, p);
    gamma = omega1 / (lambda * omega_p + omega_q);
    inv_omega2 = 1 / omega1;

    memcpy(r, beta, sizeof(double)*n);
    dscal_(&n, &negLambda, r, &inc);

    daxpy_(&n, &gamma, p, &inc, beta, &inc);
    dscal_(&active, &gamma, q, &inc);

    cglsFun2_nrOne(active, J, C, set2,
	     n0, n, q, o, z, r);

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
        << tictoc.time() << " seconds." << endl;
  }
  delete[] z;
  delete[] q;
  delete[] r;
  delete[] p;
  delete[] set2;
  return optimality;
}

int L2_SVM_MFN_nrOne(const AlgIn& data, struct options* Options,
		     struct vector_double* Weights,
		     struct vector_double* Outputs) {
  /* Disassemble the structures */
  timer tictoc;
  tictoc.restart();
  double** set = data.vals; // array of values
  const double* Y = data.Y; // labels
  const double* C = data.C; /* cost associated with each example */
  int n = Weights->d;
  int m = data.m; // number of data points
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
  int i = 0;
  vector_int* ActiveSubset = new vector_int[1];
  ActiveSubset->vec = new int[m];
  ActiveSubset->d = m;

  // double* set2 = new double[active * n];
  // double rowStart = 0;
  // initialize
  F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
  int active = 0;
  int inactive = m - 1; // l-1
  for (i = 0; i < m; i++) {
    diff = 1 - Y[i] * o[i];
    if (diff > 0) {
      ActiveSubset->vec[active] = i;
      active++;
      F += 0.5 * C[i] * diff * diff;
    } else {
      ActiveSubset->vec[inactive] = i;
      inactive--;
    }
    // // copy values over to set2
    // rowStart = i * n;
    // memcpy(set2 + rowStart, set[i], sizeof(double)*n);
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
    opt = CGLS(data,
               lambda,
               cgitermax,
               epsilon,
               ActiveSubset,
               Weights_bar,
               Outputs_bar);
    for (i = active; i < m; i++) {
      ii = ActiveSubset->vec[i];
      o_bar[ii] = ddot_(&n0, set[ii], &inc, w_bar, &inc) + w_bar[n - 1];
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
    delta = line_search_nrOne(w, w_bar, lambda, o, o_bar, Y, C, n, m);
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
        ActiveSubset->vec[active] = i;
        active++;
        F += 0.5 * C[i] * diff * diff;
      } else {
        ActiveSubset->vec[inactive] = i;
        inactive--;
      }
    }

    ActiveSubset->d = active;
    if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
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
  return 0;
}

static double l2SvmMfnFun1(double* diffs, int m, const double* Y, double* o, const double* C){
  int i = 0;
  double F = 0.0;

  #pragma omp parallel for private(i) reduction(+:F) schedule(guided)
  for (i = 0; i < m; i++) {
    double diff = 1 - Y[i] * o[i];
    diffs[i] = diff;
    if (diff > 0) {
      F += 0.5 * C[i] * diff * diff;
    }
  }

  return(F);
}

static double l2SvmMfnFun2(double* diffs, int m, const double* Y, 
    double* o, double* o_bar, double delta, const double* C){
  int i = 0;
  double F = 0.0;

  #pragma omp parallel for private(i) reduction(+:F) schedule(guided)
  for (i = 0; i < m; i++) {
  o[i] += delta * (o_bar[i] - o[i]);
    double diff = 1 - Y[i] * o[i];
    diffs[i] = diff;
    if (diff > 0) {
      F += 0.5 * C[i] * diff * diff;
    }
  }

  return(F);
}

double cglsFun1(int active, int* J, const double* C, double* set, int n0,
		int n, double* q, double* p){
  double omega_q = 0.0;
  int inc = 1;
  int i = 0;

  #pragma omp parallel for private(i) reduction(+:omega_q) schedule(guided)
  for (i = 0; i < active; i++) {
  q[i] = ddot_(&n, set + i * n, &inc, p, &inc); // + p[n-1];
    omega_q += C[J[i]] * (q[i]) * (q[i]);
  }

  return(omega_q);
}

static void cglsFun2(int active, int* J, const double* C, double* set, 
		     int n0, int n, double* q, double* o, double* z, double* r, 
		     Reduce_Vectors *reduce_vectors){
  int i;
  int inc = 1;
  int ind = 0;
  double temp = 0.0;

  reduce_vectors->init();
#pragma omp parallel for private(i) schedule(guided)
  for (i = 0; i < active; i++) {
    o[J[i]] += q[i];
    z[i] -= C[J[i]] * q[i];
    reduce_vectors->sum_scale_x(z[i], set + i * n);
  }
  reduce_vectors->reduce_sum_cgls(r);
}

static void cglsFun0(int active, int* J, const double* C, 
		     const double* Y, double* set, 
		     int n0, int n, double* o, double* z, double* r, 
		     Reduce_Vectors *reduce_vectors){
  int i;
  int ii = 0;
  int inc = 1;

  reduce_vectors->init();
  #pragma omp parallel for private(i) schedule(guided)
  for (i = 0; i < active; i++) {
    z[i] = C[J[i]] * (Y[J[i]] - o[J[i]]);
    reduce_vectors->sum_scale_x(z[i], set + i * n);
  }
  reduce_vectors->reduce_sum_cgls(r);
}

static int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
		const double epsilon, const struct vector_int* Subset,
		struct vector_double* Weights, struct vector_double* Outputs, 
		Reduce_Vectors *reduce_vectors) {
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
  const double* C = data.C;
  int n = data.n;
  //  int m  = pSet->size();
  double* beta = Weights->vec;
  double* o = Outputs->vec;
  // initialize z
  double* z = new double[active];
  double* q = new double[active];
  int ii = 0;
  int i, j;
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
    rowStart = i * n;
    memcpy(set2 + rowStart, set[J[i]], sizeof(double)*n0);
    set2[rowStart + n0] = 1.0;
  }

  cglsFun0(active, J, C, Y, set2, n0, n, o, z, r, reduce_vectors);

  double* p = new double[n];
  daxpy_(&n, &negLambda, beta, &inc, r, &inc);
  memcpy(p, r, sizeof(double)*n);
  double omega1 = ddot_(&n, r, &inc, r, &inc);
  double omega_p = omega1;
  double omega_q = 0.0;
  // double oq = 0.0;
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
    // omega_q = 0.0;
    double t = 0.0;
    omega_q = cglsFun1(active, J, C, set2, n0, n, q, p);

    gamma = omega1 / (lambda * omega_p + omega_q);
    inv_omega2 = 1 / omega1;

    memcpy(r, beta, sizeof(double)*n);
    dscal_(&n, &negLambda, r, &inc);

    daxpy_(&n, &gamma, p, &inc, beta, &inc);
    dscal_(&active, &gamma, q, &inc);

    cglsFun2(active, J, C, set2,
	     n0, n, q, o, z, r, reduce_vectors);

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
               struct vector_double* Outputs) {
  /* Disassemble the structures */
  timer tictoc;
  tictoc.restart();
  double** set = data.vals; // array of values
  const double* Y = data.Y; // labels
  const double* C = data.C; /* cost associated with each example */
  int n = Weights->d;
  int m = data.m; // number of data points
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
  int i = 0;
  vector_int* ActiveSubset = new vector_int[1];
  ActiveSubset->vec = new int[m];
  ActiveSubset->d = m;

  int nr_thread = THREADS;
  omp_set_num_threads(nr_thread);

  // Need accumulators for OMP
  Reduce_Vectors *reduce_vectors = new Reduce_Vectors(n);

  // initialize
  F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);

  int active = 0;
  int inactive = m - 1; // l-1
#ifdef FUN2
  double* diffs = new double[m];
  F = F + l2SvmMfnFun1(diffs, m, Y, o, C);
  for (i = 0; i < m; i++) {
    if (diffs[i] > 0) {
      ActiveSubset->vec[active] = i;
      active++;
    } else {
      ActiveSubset->vec[inactive] = i;
      inactive--;
    }
  }
#else
  for (i = 0; i < m; i++) {
    diff = 1 - Y[i] * o[i];
    if (diff > 0) {
      ActiveSubset->vec[active] = i;
      active++;
      F += 0.5 * C[i] * diff * diff;
    } else {
      ActiveSubset->vec[inactive] = i;
      inactive--;
    }
  }
#endif

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
    opt = CGLS(data,
               lambda,
               cgitermax,
               epsilon,
               ActiveSubset,
               Weights_bar,
               Outputs_bar, 
	       reduce_vectors);
    for (i = active; i < m; i++) {
      ii = ActiveSubset->vec[i];
      o_bar[ii] = ddot_(&n0, set[ii], &inc, w_bar, &inc) + w_bar[n - 1];
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
        delete[] ActiveSubset->vec;
        delete[] ActiveSubset;
        delete[] o_bar;
        delete[] w_bar;
#ifdef FUN2
	delete[] diffs;
#endif
        // delete[] o_barD;
        // delete[] w_barD;
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
    // for (int i = n; i--;) {
    //   w_barD[i] = w_bar[i] - w[i];
    // }
    // for (int i = m; i--;) {
    //   o_barD[i] = o_bar[i] - o[i];
    // }

    delta = line_search(w, w_bar, lambda, o, o_bar, Y, C, n, m); //, w_barD, o_barD);
    F_old = F;
    // daxpy_(&n, &delta, w_barD, &inc, w, &inc);
    double delta2 = 1-delta;
    dscal_(&n, &delta2, w, &inc);
    daxpy_(&n, &delta, w_bar, &inc, w, &inc);
    F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);

    active = 0;
    inactive = m - 1;

#ifdef FUN2
    F = F + l2SvmMfnFun2(diffs, m, Y, 
    			 o, o_bar, delta, C);
    for (int i = 0; i < m; i++) {
      if (diffs[i] > 0) {
        ActiveSubset->vec[active] = i;
        active++;
      } else {
        ActiveSubset->vec[inactive] = i;
        inactive--;
      }
    }
#else
    for (int i = 0; i < m; i++) {
      o[i] += delta * (o_bar[i] - o[i]);
      diff = 1 - Y[i] * o[i];
      if (diff > 0) {
        ActiveSubset->vec[active] = i;
        active++;
        F += 0.5 * C[i] * diff * diff;
      } else {
        ActiveSubset->vec[inactive] = i;
        inactive--;
      }
    }
#endif
    ActiveSubset->d = active;
    if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
      return 2;
    }
  }
  delete[] ActiveSubset->vec;
  delete[] ActiveSubset;
  delete[] o_bar;
  delete[] w_bar;
#ifdef FUN2
  delete[] diffs;
#endif
  // delete[] o_barD;
  // delete[] w_barD;
  delete[] Weights_bar;
  delete[] Outputs_bar;
  delete reduce_vectors;
  tictoc.stop();
  return 0;
}

static void lsFun1(int l, double& L, double& R, 
		   const double* Y, double* o, 
		   double* o_bar,
		   const double* C){
  int i = 0;
  double Lt = 0.0;
  double Rt = 0.0;

  #pragma omp parallel for private(i) reduction(+:Lt, Rt) schedule(guided)
  for (i = 0; i < l; i++) {
    if (Y[i] * o[i] < 1) {
      double diff = C[i] * (o_bar[i] - o[i]);
      Lt += (o[i] - Y[i]) * diff;
      Rt += (o_bar[i] - Y[i]) * diff;
    }
  }
  L += Lt;
  R += Rt;
}

double line_search(double* w, double* w_bar, double lambda, double* o,
                   double* o_bar, const double* Y, const double* C, int d, /* data dimensionality -- 'n' */
                   int l){ // , double* w_barD, double* o_barD
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

#ifdef LFUN2
  lsFun1(l, L, R, Y, o, o_bar, C);
#endif
  // for (i = 0; i < l; i++) {
  //   if (Y[i] * o[i] < 1) {
  //     diff = C[i] * (o_bar[i] - o[i]);
  //     L += (o[i] - Y[i]) * diff;
  //     R += (o_bar[i] - Y[i]) * diff;
  //   }
  // }
  Delta* deltas = new Delta[l];
  int p = 0;
  for (i = 0; i < l; i++) {
    diff = Y[i] * (o_bar[i] - o[i]);
    if (Y[i] * o[i] < 1) {
#ifndef LFUN2
      d2 = C[i] * (o_bar[i] - o[i]);
      L += (o[i] - Y[i]) * d2;
      R += (o_bar[i] - Y[i]) * d2;
#endif
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
    diff = (deltas[i].s) * C[ii] * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  delete[] deltas;
  return (-L / (R - L));
}
 
double line_search_nrOne(double* w, double* w_bar, double lambda, double* o,
			 double* o_bar, const double* Y, const double* C, int d, /* data dimensionality -- 'n' */
			 int l){ // , double* w_barD, double* o_barD
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
      d2 = C[i] * (o_bar[i] - o[i]);
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
    diff = (deltas[i].s) * C[ii] * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  delete[] deltas;
  return (-L / (R - L));
}

/********************** TRON, heavily modified from liblinear implementation, 
                        optimized for single and multi-threaded use within Percolator
Author: John T. Halloran
Affiliation: UC Davis
Date: May 2017
********************/

static void default_print(const char *buf)
{
	fputs(buf,stdout);
	fflush(stdout);
}

void troninfo(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	default_print(buf);
}

static void axpy(const int n, const double a, const double* x, double *y)
{
  for(int i = 0; i < n; i++){
    y[i] += a * x[i];
  }
}

static double norm_inf(int n, double *x)
{
	double dmax = fabs(x[0]);
	for (int i=1; i<n; i++)
		if (fabs(x[i]) >= dmax)
			dmax = fabs(x[i]);
	return(dmax);
}


double fun_omp(double *w, int w_size, int l,
	       const double* y, double* z, const double* C, double* X)
{
  unsigned int i;
  int inc = 1;
  double f=0;
// #ifdef BLASC
// f = cblas_ddot(w_size, w, 1, w, 1);
// #else
f = ddot_(&w_size, w, &inc, w, &inc) / 2.0;
// #endif
  f /= 2.0;

  // #pragma omp parallel for private(i) reduction(+:f) schedule(static)
  #pragma omp parallel for private(i) reduction(+:f) schedule(guided)
  for(i=0;i<l;i++){
  z[i]=y[i]*ddot_(&w_size, w, &inc, X + i * w_size, &inc); // + w[w_size - 1]);
// #ifdef BLASC
//       z[i]=y[i]*(cblas_ddot(w_size, w, 1, X[i], 1) + w[w_size - 1]);
// #else
  
      // #endif

      double d = 1-z[i];
      if (d > 0)
        f += C[i]*d*d;
    }

  return(f);
}

static void subXTv(double *v, double *XTv, int w_size,
    double* X, const int* Id, int sizeI, Reduce_Vectors *reduce_vectors)
{
  int i;
  int inc = 1;
  double a = 0;

  reduce_vectors->init();

#pragma omp parallel for private(i) schedule(guided)
  for(i=0;i<sizeI;i++)
    reduce_vectors->sum_scale_x(v[i], X + Id[i] * w_size);
  
  reduce_vectors->reduce_sum(XTv);
}

static int grad(double *w, double *g, int w_size, int l,
		const double* y, double* z, const double* C, double* X,
		int* Id, Reduce_Vectors *reduce_vectors)
{
  int i;
  int sizeI = 0;
        
  for (i=0;i<l;i++){
    if (z[i] < 1) {
      z[sizeI] = C[i]*y[i]*(z[i]-1);
      Id[sizeI] = i;
      sizeI++;
    }
  }

  subXTv(z, g, w_size, X, Id, sizeI, reduce_vectors);

  for(i=0;i<w_size;i++)
    g[i] = w[i] + 2*g[i];
  return(sizeI);
}

static void Hv(double *s, double *Hs, int w_size,
	       const double* C, double* X, const int *Id, int sizeI,
	       Reduce_Vectors *reduce_vectors)
{
  int i;
  int inc = 1;

  reduce_vectors->init();
#pragma omp parallel for private(i) schedule(guided)
  for(i=0;i<sizeI;i++)
    {
      double* xi = X + Id[i] * w_size;
// #ifdef BLASC
//       double xTs = C[Id[i]]*(cblas_ddot(w_size, s, inc, xi, inc) + s[w_size - 1]);
// #else
      double xTs = C[Id[i]]*ddot_(&w_size, s, &inc, xi, &inc); // + s[w_size - 1]);
      // #endif
      reduce_vectors->sum_scale_x(xTs, xi);
    }
  reduce_vectors->reduce_sum(Hs);
  for(i=0;i<w_size;i++)
    Hs[i] = s[i] + 2*Hs[i];
}

static int trcg(double delta, double *g, double *s, double *r, bool *reach_boundary, 
		int n, const double* C, double* X, const int* Id, int sizeI, 
		const int cgitermax, 
		Reduce_Vectors *reduce_vectors,
		const double eps_cg = 0.1)
{
  int i, inc = 1;
  double one = 1;
  double *d = new double[n];
  double *Hd = new double[n];
  double rTr, rnewTrnew, alpha, beta, cgtol;

  *reach_boundary = false;
  for (i=0; i<n; i++)
    {
      s[i] = 0;
      r[i] = -g[i];
      d[i] = r[i];
    }
// #ifdef BLASC
//   cgtol = eps_cg*cblas_dnrm2(n, g, inc);
// #else
  cgtol = eps_cg*dnrm2_(&n, g, &inc);
  // #endif

  int cg_iter = 0;
// #ifdef BLASC
//   rTr = cblas_ddot(n, r, inc, r, inc);
// #else
  rTr = ddot_(&n, r, &inc, r, &inc);
  // #endif
  // while (cg_iter < cgitermax)
  while (1)
    {
// #ifdef BLASC
//       if (cblas_dnrm2(n, r, inc) <= cgtol) break;
// #else
      if (dnrm2_(&n, r, &inc) <= cgtol) break;
      // #endif
      cg_iter++;
      Hv(d, Hd, n,
         C, X, Id, sizeI,
	 reduce_vectors);

// #ifdef BLASC
//       alpha = rTr/cblas_ddot(n, d, inc, Hd, inc);
//       cblas_daxpy(n, alpha, d, inc, s, inc);
//       if (cblas_dnrm2(n, s, inc) > delta)
// #else
	alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
	daxpy_(&n, &alpha, d, &inc, s, &inc);
      if (dnrm2_(&n, s, &inc) > delta)
	// #endif
        {
          // troninfo("cg reaches trust region boundary\n");
          *reach_boundary = true;
          alpha = -alpha;
// #ifdef BLASC
//           cblas_daxpy(n, alpha, d, inc, s, inc);

//           double std = cblas_ddot(n, s, inc, d, inc);
//           double sts = cblas_ddot(n, s, inc, s, inc);
//           double dtd = cblas_ddot(n, d, inc, d, inc);
// #else
	  daxpy_(&n, &alpha, d, &inc, s, &inc);

	  double std = ddot_(&n, s, &inc, d, &inc);
	  double sts = ddot_(&n, s, &inc, s, &inc);
	  double dtd = ddot_(&n, d, &inc, d, &inc);
	  // #endif
          double dsq = delta*delta;
          double rad = sqrt(std*std + dtd*(dsq-sts));
          if (std >= 0)
            alpha = (dsq - sts)/(std + rad);
          else
            alpha = (rad - std)/dtd;
// #ifdef BLASC
//           cblas_daxpy(n, alpha, d, inc, s, inc);
//           alpha = -alpha;
//           cblas_daxpy(n, alpha, Hd, inc, r, inc);
// #else
	  daxpy_(&n, &alpha, d, &inc, s, &inc);
	  alpha = -alpha;
	  daxpy_(&n, &alpha, Hd, &inc, r, &inc);
	  //#endif
          break;
        }
      alpha = -alpha;
// #ifdef BLASC
//       cblas_daxpy(n, alpha, Hd, inc, r, inc);
//       rnewTrnew = cblas_ddot(n, r, inc, r, inc);
//       beta = rnewTrnew/rTr;
//       cblas_dscal(n, beta, d, inc);
//       cblas_daxpy(n, one, r, inc, d, inc);
// #else
      daxpy_(&n, &alpha, Hd, &inc, r, &inc);
      rnewTrnew = ddot_(&n, r, &inc, r, &inc);
      beta = rnewTrnew/rTr;
      dscal_(&n, &beta, d, &inc);
      daxpy_(&n, &one, r, &inc, d, &inc);
      // #endif
      rTr = rnewTrnew;
    }

  delete[] d;
  delete[] Hd;

  return(cg_iter);
}

int tron(const AlgIn& data, struct options* Options,
         struct vector_double* Weights,
         struct vector_double* Outputs)
{
  timer tictoc;
  tictoc.restart();

  int nr_thread = THREADS;
  omp_set_num_threads(nr_thread);

  // Parameters for updating the iterates.
  double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;

  // Parameters for updating the trust region size delta.
  double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

  double eps = 0.01;
  double eps_cg = 0.1;
  int max_iter = Options->mfnitermax;
  // int max_iter = 10;
  int cgitermax = SMALL_CGITERMAX;
  double* w = Weights->vec; // weights to be learned
  const double* Y = data.Y; // labels
  const double* C = data.C; /* cost associated with each example */
  int n = Weights->d;
  int l = data.m;
  int* Id = new int[l];
  int i, cg_iter;
  double delta, snorm, one=1.0;
  double alpha, f, fnew, prered, actred, gs;
  int search = 1, iter = 1, inc = 1;
  double* s = new double[n];
  double* r = new double[n];
  double* g = new double[n];
  double* z = new double[l];
  int sizeI;
  int ini = 0;

  double** X0 = data.vals; // array of values

  double* X = new double[l * n];
  unsigned int ind = 0;
  unsigned int numZeros = 0;

  for(i = 0; i < l; i++){
    memcpy(X+ind, X0[i], sizeof(double)*n);
    X[ind+n-1] = 1.0;
    ind += n;
  }

  // Need accumulators for OMP
  Reduce_Vectors *reduce_vectors = new Reduce_Vectors(n);

  // /////////////
  // double* o = Outputs->vec; // predictions of w on x
  // if the above is needed, set: double* o = z; and update the arrays accordingly

  f = fun_omp(w, n, l,
	      Y, z, C, X);
  sizeI = grad(w, g, n, l,
               Y, z, C, X, Id, reduce_vectors);

// #ifdef BLASC
//   delta = cblas_dnrm2(n, g, inc);
// #else
  delta = dnrm2_(&n, g, &inc);
  //#endif
  double gnorm = delta;

  if (gnorm < eps)
    search = 0;

  iter = 1;

  double *w_new = new double[n];
  bool reach_boundary;
  while (iter <= max_iter && search)
    {
      cg_iter = trcg(delta, g, s, r, &reach_boundary, n, 
		     C, X, Id, sizeI, cgitermax, 
		     reduce_vectors,
		     eps_cg);

      memcpy(w_new, w, sizeof(double)*n);
// #ifdef BLASC
//       cblas_daxpy(n, one, s, inc, w_new, inc);

//       gs = cblas_ddot(n, g, inc, s, inc);
//       prered = -0.5*(gs-cblas_ddot(n, s, inc, r, inc));
// #else
      daxpy_(&n, &one, s, &inc, w_new, &inc);

      gs = ddot_(&n, g, &inc, s, &inc);
      prered = -0.5*(gs-ddot_(&n, s, &inc, r, &inc));
      // #endif

      fnew = fun_omp(w_new, n, l,
		     Y, z, C, X);
      // Compute the actual reduction.
      actred = f - fnew;

      // On the first iteration, adjust the initial step bound.
// #ifdef BLASC
//       snorm = cblas_dnrm2(n, s, inc);
// #else
      snorm = dnrm2_(&n, s, &inc);
      // #endif
      if (iter == 1)
        delta = min(delta, snorm);

      // Compute prediction alpha*snorm of the step.
      if (fnew - f - gs <= 0)
        alpha = sigma3;
      else
        alpha = max(sigma1, -0.5*(gs/(fnew - f - gs)));

      // Update the trust region bound according to the ratio of actual to predicted reduction.
      if (actred < eta0*prered)
        delta = min(max(alpha, sigma1)*snorm, sigma2*delta);
      else if (actred < eta1*prered)
        delta = max(sigma1*delta, min(alpha*snorm, sigma2*delta));
      else if (actred < eta2*prered)
        delta = max(sigma1*delta, min(alpha*snorm, sigma3*delta));
      else
        {
          if (reach_boundary)
            delta = sigma3*delta;
          else
            delta = max(delta, min(alpha*snorm, sigma3*delta));
        }

      // troninfo("iter %2d act %5.3e pre %5.3e delta %5.3e f %5.3e |g| %5.3e CG %3d\n", iter, actred, prered, delta, f, gnorm, cg_iter);

      if (actred > eta0*prered)
        {
          iter++;
          memcpy(w, w_new, sizeof(double)*n);
          f = fnew;
          sizeI = grad(w, g, n, l,
                       Y, z, C, X, Id, reduce_vectors);
// #ifdef BLASC
//           gnorm = cblas_dnrm2(n, g, inc);
// #else
	  gnorm = dnrm2_(&n, g, &inc);
	  // #endif

	  if (gnorm <= eps)
	    break;
        }
      if (f < -1.0e+32)
        {
          // troninfo("WARNING: f < -1.0e+32\n");
          break;
        }
      if (prered <= 0)
        {
          // troninfo("WARNING: prered <= 0\n");
          break;
        }
      if (fabs(actred) <= 1.0e-12*fabs(f) &&
          fabs(prered) <= 1.0e-12*fabs(f))
        {
          // troninfo("WARNING: actred and prered too small\n");
          break;
        }
    }

  delete[] g;
  delete[] r;
  delete[] w_new;
  delete[] s;
  delete[] Id;
  delete [] z; 
  delete [] X;
  delete reduce_vectors;

  tictoc.stop();
  return 0;
}

///////////////////////////////////////
int tron_nrOne(const AlgIn& data, struct options* Options,
	       struct vector_double* Weights,
	       struct vector_double* Outputs)
{
  timer tictoc;
  tictoc.restart();

  // openblas_set_num_threads(1);
  // Parameters for updating the iterates.
  double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;

  // Parameters for updating the trust region size delta.
  double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

  double eps = 0.01;
  double eps_cg = 0.1;
  int max_iter = Options->mfnitermax;
  // int max_iter = 10;
  int cgitermax = SMALL_CGITERMAX;
  double* w = Weights->vec; // weights to be learned
  const double* Y = data.Y; // labels
  const double* C = data.C; /* cost associated with each example */
  int n = Weights->d;
  int l = data.m;
  int* Id = new int[l];
  int i, cg_iter;
  double delta, snorm, one=1.0;
  double alpha, f, fnew, prered, actred, gs;
  int search = 1, iter = 1, inc = 1;
  double* s = new double[n];
  double* r = new double[n];
  double* g = new double[n];
  double* z = new double[l];
  double* cProdY = new double[l];
  int sizeI;
  int ini = 0;

  // Load feature matrix into contiguous memory as flat array
  ////////////// features
  double** X0 = data.vals; // array of values

  double* X = new double[l * n];
  unsigned int ind = 0;
  unsigned int numZeros = 0;

  for(i = 0; i < l; i++){
    memcpy(X+ind, X0[i], sizeof(double)*n);
    X[ind+n-1] = 1.0;
    ind += n;
  }
  // cout << numZeros << " zero entries of " << ind << " total features\n";

  for(i = 0; i < l; i++){
    cProdY[i] = C[i] * Y[i];
  }

  f = fun(w, n, l,
          Y, z, C, X);
  sizeI = grad(w, g, n, l,
	       cProdY, z, X, Id);
// #ifdef BLASC
//   delta = cblas_dnrm2(n, g, inc);
// #else
  delta = dnrm2_(&n, g, &inc);
  // #endif
  double gnorm = delta;

  if (gnorm < eps)
    search = 0;

  iter = 1;

  double *w_new = new double[n];
  bool reach_boundary;
  while (iter <= max_iter && search)
    {
      cg_iter = trcg(delta, g, s, r, &reach_boundary, n, 
		     C, X, Id, sizeI, cgitermax, eps_cg);

      memcpy(w_new, w, sizeof(double)*n);
// #ifdef BLASC
//       cblas_daxpy(n, one, s, inc, w_new, inc);

//       gs = cblas_ddot(n, g, inc, s, inc);
//       prered = -0.5*(gs-cblas_ddot(n, s, inc, r, inc));
// #else
      daxpy_(&n, &one, s, &inc, w_new, &inc);
      gs = ddot_(&n, g, &inc, s, &inc);
      prered = -0.5*(gs-ddot_(&n, s, &inc, r, &inc));
      // #endif

      fnew = fun(w_new, n, l,
                 Y, z, C, X);

      // Compute the actual reduction.
      actred = f - fnew;

      // On the first iteration, adjust the initial step bound.
// #ifdef BLASC
//       snorm = cblas_dnrm2(n, s, inc);
// #else
      snorm = dnrm2_(&n, s, &inc);
      // #endif
      if (iter == 1)
        delta = min(delta, snorm);

      // Compute prediction alpha*snorm of the step.
      if (fnew - f - gs <= 0)
        alpha = sigma3;
      else
        alpha = max(sigma1, -0.5*(gs/(fnew - f - gs)));

      // Update the trust region bound according to the ratio of actual to predicted reduction.
      if (actred < eta0*prered)
        delta = min(max(alpha, sigma1)*snorm, sigma2*delta);
      else if (actred < eta1*prered)
        delta = max(sigma1*delta, min(alpha*snorm, sigma2*delta));
      else if (actred < eta2*prered)
        delta = max(sigma1*delta, min(alpha*snorm, sigma3*delta));
      else
        {
          if (reach_boundary)
            delta = sigma3*delta;
          else
            delta = max(delta, min(alpha*snorm, sigma3*delta));
        }

      if (actred > eta0*prered)
        {
          iter++;
          memcpy(w, w_new, sizeof(double)*n);
          f = fnew;
          sizeI = grad(w, g, n, l,
		       cProdY, z, X, Id);
// #ifdef BLASC
//           gnorm = cblas_dnrm2(n, g, inc);
// #else
          gnorm = dnrm2_(&n, g, &inc);
	  // #endif

	  if (gnorm <= eps)
	    break;
        }
      if (f < -1.0e+32)
        {
          // troninfo("WARNING: f < -1.0e+32\n");
          break;
        }
      if (prered <= 0)
        {
          // troninfo("WARNING: prered <= 0\n");
          break;
        }
      if (fabs(actred) <= 1.0e-12*fabs(f) &&
          fabs(prered) <= 1.0e-12*fabs(f))
        {
          break;
        }
    }

  delete[] g;
  delete[] r;
  delete[] w_new;
  delete[] s;
  delete[] Id;
  delete[] z; 
  delete[] cProdY;
  delete [] X;
  tictoc.stop();

  return 0;
}

int trcg(double delta, double *g, double *s, double *r, bool *reach_boundary, 
	 int n, const double* C, double* X, const int* Id, int sizeI, 
	 const int cgitermax, const double eps_cg = 0.1)
{
  int i, inc = 1;
  double one = 1;
  double *d = new double[n];
  double *Hd = new double[n];
  double *Xs = new double[sizeI];
  double rTr, rnewTrnew, alpha, beta, cgtol;

  *reach_boundary = false;
  for (i=0; i<n; i++)
    {
      s[i] = 0;
      r[i] = -g[i];
      d[i] = r[i];
    }
// #ifdef BLASC
//   cgtol = eps_cg*cblas_dnrm2(n, g, inc);

//   int cg_iter = 0;
//   rTr = cblas_ddot(n, r, inc, r, inc);
// #else
  cgtol = eps_cg*dnrm2_(&n, g, &inc);

  int cg_iter = 0;
  rTr = ddot_(&n, r, &inc, r, &inc);

  //#endif

  // while (cg_iter < cgitermax)
  while (1)
    {
// #ifdef BLASC
//       if (cblas_dnrm2(n, r, inc) <= cgtol)
// #else
      if (dnrm2_(&n, r, &inc) <= cgtol)
	// #endif
        break;
      cg_iter++;
      Hv(d, Hd, n,
         C, X, Id, sizeI);

// #ifdef BLASC      
//       alpha = rTr/cblas_ddot(n, d, inc, Hd, inc);
//       cblas_daxpy(n, alpha, d, inc, s, inc);
//       if (cblas_dnrm2(n, s, inc) > delta)
//         {
//           // troninfo("cg reaches trust region boundary\n");
//           *reach_boundary = true;
//           alpha = -alpha;
//           cblas_daxpy(n, alpha, d, inc, s, inc);

//           double std = cblas_ddot(n, s, inc, d, inc);
//           double sts = cblas_ddot(n, s, inc, s, inc);
//           double dtd = cblas_ddot(n, d, inc, d, inc);
//           double dsq = delta*delta;
//           double rad = sqrt(std*std + dtd*(dsq-sts));
//           if (std >= 0)
//             alpha = (dsq - sts)/(std + rad);
//           else
//             alpha = (rad - std)/dtd;
//           cblas_daxpy(n, alpha, d, inc, s, inc);
//           alpha = -alpha;
//           cblas_daxpy(n, alpha, Hd, inc, r, inc);
//           break;
//         }
//       alpha = -alpha;
//       cblas_daxpy(n, alpha, Hd, inc, r, inc);
//       rnewTrnew = cblas_ddot(n, r, inc, r, inc);
//       beta = rnewTrnew/rTr;
//       cblas_dscal(n, beta, d, inc);
//       cblas_daxpy(n, one, r, inc, d, inc);
//       rTr = rnewTrnew;
//     }
// #else
      alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
      daxpy_(&n, &alpha, d, &inc, s, &inc);
      if (dnrm2_(&n, s, &inc) > delta)
        {
          // troninfo("cg reaches trust region boundary\n");
          *reach_boundary = true;
          alpha = -alpha;
          daxpy_(&n, &alpha, d, &inc, s, &inc);

          double std = ddot_(&n, s, &inc, d, &inc);
          double sts = ddot_(&n, s, &inc, s, &inc);
          double dtd = ddot_(&n, d, &inc, d, &inc);
          double dsq = delta*delta;
          double rad = sqrt(std*std + dtd*(dsq-sts));
          if (std >= 0)
            alpha = (dsq - sts)/(std + rad);
          else
            alpha = (rad - std)/dtd;
          daxpy_(&n, &alpha, d, &inc, s, &inc);
          alpha = -alpha;
          daxpy_(&n, &alpha, Hd, &inc, r, &inc);
          break;
        }
      alpha = -alpha;
      daxpy_(&n, &alpha, Hd, &inc, r, &inc);
      rnewTrnew = ddot_(&n, r, &inc, r, &inc);
      beta = rnewTrnew/rTr;
      dscal_(&n, &beta, d, &inc);
      daxpy_(&n, &one, r, &inc, d, &inc);
      rTr = rnewTrnew;
    }
  // #endif

  delete[] d;
  delete[] Hd;
  return(cg_iter);
}

double fun(double *w, int w_size, int l,
	   const double* y, double* z, const double* C, double* X)
{
  int i;
  double f=0;
  double d = 0;
  int inc = 1;
  
// #ifdef BLASC
//   f = cblas_ddot(w_size, w, 1, w, 1);
// #else
  f = ddot_(&w_size, w, &inc, w, &inc) / 2.0;
// #endif

// #ifdef BLASC
//   cblas_dgemv(CblasRowMajor, CblasNoTrans,
//   	      l, w_size, 1.0, X, w_size,
//   	      w, 1, 0.0, z, 1);
// #else
  char trans = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  dgemv_(&trans, &w_size, &l,
	 &alpha, X, &w_size,
	 w, &inc, &beta, z, &inc);
  // #endif

  for(i=0;i<l;i++)
    {
      z[i]=y[i]*z[i];
      d = 1-z[i];
      if (d > 0)
        f += C[i]*d*d;
    }

  return(f);
}

int grad(double *w, double *g, int w_size, int l,
	 const double* cProdY, double* z, double* X,
	 int* Id)
{
  int i;
  int sizeI = 0;

  for (i=0;i<l;i++){
    if (z[i] < 1) {
      z[sizeI] = cProdY[i]*(z[i]-1);
      Id[sizeI] = i;
      sizeI++;
    }
  }

  subXTv(z, g, w_size, X, Id, sizeI);

  for(i=0;i<w_size;i++)
    g[i] = w[i] + 2*g[i];
  return(sizeI);
}

void Hv(double *s, double *Hs, int w_size,
        const double* C, double* X, const int *Id, int sizeI)
{
  int i;
  int inc = 1;
  unsigned int ind = 0;
  double xTs = 0;

  for(i=0;i<w_size;i++)
    Hs[i]=0;
  for(i=0;i<sizeI;i++)
    {

// #ifdef BLASC
//       xTs = C[Id[i]] * (cblas_ddot(w_size, s, inc, X + Id[i] * w_size, inc) + s[w_size - 1]);
//       cblas_daxpy(w_size, xTs, X + Id[i] * w_size, inc, Hs, inc);
// #else
      xTs = C[Id[i]] * ddot_(&w_size, s, &inc, X + Id[i] * w_size, &inc); // + s[w_size - 1]);
      daxpy_(&w_size, &xTs, X + Id[i] * w_size, &inc, Hs, &inc);
      // #endif
    }
  for(i=0;i<w_size;i++)
    Hs[i] = s[i] + 2*Hs[i];
}

void subXTv(double *v, double *XTv, int w_size,
	    double* X, const int* Id, int sizeI)
{
  int i;
  int inc = 1;

  for(i=0;i<w_size;i++)
    XTv[i]=0;
  for(i=0;i<sizeI;i++)
    {
// #ifdef BLASC
//       cblas_daxpy(w_size, v[i], X + Id[i] * w_size, inc, XTv, inc);
// #else
      daxpy_(&w_size, &(v[i]), X + Id[i] * w_size, &inc, XTv, &inc);
      // #endif
    }

}
