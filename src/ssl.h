/*******************************************************************************
 * SVMlin
 * Copyright (c) 2006 Vikas Sindhwani at the University of Chicago.
 * Adapted to Percolator by Lukas Käll at the University of Washington
 *
 *******************************************************************************/
#ifndef _svmlin_H
#define _svmlin_H
#include <vector>
#include <ctime>

using namespace std;

/* OPTIMIZATION CONSTANTS */
#define CGITERMAX 10000 /* maximum number of CGLS iterations */
#define SMALL_CGITERMAX 10 /* for heuristic 1 in reference [2] */
#define EPSILON   1e-7 /* most tolerances are set to this value */
#define BIG_EPSILON 0.01 /* for heuristic 2 in reference [2] */
#define RELATIVE_STOP_EPS 1e-9 /* for L2-SVM-MFN relative stopping criterion */
#define MFNITERMAX 50 /* maximum number of MFN iterations */

#define VERBOSE_CGLS 0

class AlgIn {
  public:
    AlgIn(const int size, const int numFeat);
    virtual ~AlgIn();
    int m; /* number of examples */
    int n; /* number of features */
    int positives;
    int negatives;
    const double** vals;
    double* Y; /* labels */
    double* C; /* cost associated with each example */
    void setCost(double pos, double neg) {
      int ix = 0;
      for (; ix < negatives; ++ix) {
        C[ix] = neg;
      }
      for (; ix < negatives + positives; ++ix) {
        C[ix] = pos;
      }
    }
};

/* Data: Input examples are stored in sparse (Compressed Row Storage) format */
struct data {
    int m; /* number of examples */
    int l; /* number of labeled examples */
    int u; /* number of unlabeled examples l+u = m */
    int n; /* number of features */
    int nz; /* number of non-zeros */
    double* val; /* data values (nz elements) [CRS format] */
    int* rowptr; /* n+1 vector [CRS format] */
    int* colind; /* nz elements [CRS format] */
    double* Y; /* labels */
    double* C; /* cost associated with each example */
};

struct vector_double { /* defines a vector of doubles */
    int d; /* number of elements */
    double* vec; /* ptr to vector elements*/
};

struct vector_int { /* defines a vector of ints for index subsets */
    int d; /* number of elements */
    int* vec; /* ptr to vector elements */
};

struct options {
    /* user options */
    double lambda; /* regularization parameter */
    double lambda_u; /* regularization parameter over unlabeled examples */
    double epsilon; /* all tolerances */
    int cgitermax; /* max iterations for CGLS */
    int mfnitermax; /* max iterations for L2_SVM_MFN */

};

class timer { /* to output run time */
  protected:
    double start, finish;
  public:
    vector<double> times;
    void record() {
      times.push_back(time());
    }
    void reset_vectors() {
      times.erase(times.begin(), times.end());
    }
    void restart() {
      start = clock();
    }
    void stop() {
      finish = clock();
    }
    double time() const {
      return ((double)(finish - start)) / CLOCKS_PER_SEC;
    }
};
class Delta { /* used in line search */
  public:
    Delta() {
      delta = 0.0;
      index = 0;
      s = 0;
    }
    ;
    double delta;
    int index;
    int s;
};
inline bool operator<(const Delta& a, const Delta& b) {
  return (a.delta < b.delta);
}
;

void Clear(struct data* a); /* deletes a */
void Clear(struct vector_double* a); /* deletes a */
void Clear(struct vector_int* a); /* deletes a */
double norm_square(const vector_double* A); /* returns squared length of A */

/* svmlin algorithms and their subroutines */

/* Conjugate Gradient for Sparse Linear Least Squares Problems */
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_{i in Subset} Data->C[i] (Y[i]- w' x_i)^2 */
/* over a subset of examples x_i specified by vector_int Subset */
int CGLS(const AlgIn& set, const double lambda, const int cgitermax,
         const double epsilon, const struct vector_int* Subset,
         struct vector_double* Weights, struct vector_double* Outputs, double cpos, double cneg);

/* Linear Modified Finite Newton L2-SVM*/
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_i Data->C[i] max(0,1 - Y[i] w' x_i)^2 */
int L2_SVM_MFN(const AlgIn& set, struct options* Options,
               struct vector_double* Weights,
               struct vector_double* Outputs, double cpos, double cneg);
double line_search(double* w, double* w_bar, double lambda, double* o,
                   double* o_bar, const double* Y, const double* C, int d,
                   int l, double cpos, double cneg);

#endif
