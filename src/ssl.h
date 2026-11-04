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
#ifndef _svmlin_H
#define _svmlin_H
#include <cassert>
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

/* Training examples of the SVM: rows of n values, the n - 1 features
   followed by the intercept 1.0. By default the features are copied into
   one contiguous row-major matrix, so that the solver streams through
   memory instead of following pointers to scattered rows. With rowsInPlace
   only pointers to the rows are kept, which saves the copies (8 n bytes per
   example) but makes the solver slower; the rows must then hold the
   intercept, as the rows of a FeatureMemoryPool set up for this do, and
   outlive this. */
class AlgIn {
  public:
    explicit AlgIn(const int numFeat, bool rowsInPlace = false)
        : n(numFeat), rowsInPlace_(rowsInPlace) {}
    int m = 0; /* number of examples */
    int n; /* number of features, including the intercept */
    int positives = 0;
    int negatives = 0;
    std::vector<double> Y; /* labels */

    const double* row(int i) const {
      return rowsInPlace_
                 ? rowPointers_[static_cast<size_t>(i)]
                 : rows_.data() + static_cast<size_t>(i) * static_cast<size_t>(n);
    }
    bool rowsInPlace() const { return rowsInPlace_; }
    void clear() {
      m = positives = negatives = 0;
      rows_.clear();
      rowPointers_.clear();
      Y.clear();
    }
    /* room for numExamples examples, allocated but only touched when used */
    void reserve(size_t numExamples) {
      if (rowsInPlace_) {
        rowPointers_.reserve(numExamples);
      } else {
        rows_.reserve(numExamples * static_cast<size_t>(n));
      }
      Y.reserve(numExamples);
    }
    /* appends an example with the n - 1 features of row */
    void add(const double* row, double label) {
      if (rowsInPlace_) {
        assert(row[n - 1] == 1.0);
        rowPointers_.push_back(row);
      } else {
        rows_.insert(rows_.end(), row, row + (n - 1));
        rows_.push_back(1.0);
      }
      Y.push_back(label);
      ++m;
    }

  private:
    bool rowsInPlace_;
    std::vector<double> rows_; /* m x n, the copies */
    std::vector<const double*> rowPointers_; /* m, with rowsInPlace */
};

/* Asks the processor to load the n values of a row into the cache. The
   solver visits the rows of its active examples, or of the inactive ones,
   which skip through the matrix, and requests each row a few rows ahead. */
inline void prefetchRow(const double* row, int n) {
#if defined(__GNUC__) || defined(__clang__)
  const char* p = reinterpret_cast<const char*>(row);
  const int bytes = n * static_cast<int>(sizeof(double));
  for (int offset = 0; offset < bytes; offset += 64) {
    __builtin_prefetch(p + offset);
  }
  __builtin_prefetch(p + bytes - 1);
#else
  (void)row;
  (void)n;
#endif
}

struct options {
    /* user options */
    double lambda; /* regularization parameter */
    double lambda_u; /* regularization parameter over unlabeled examples */
    double epsilon; /* all tolerances */
    int cgitermax; /* max iterations for CGLS */
    int mfnitermax; /* max iterations for L2_SVM_MFN */
    /* solve the least squares problems exactly (normal equations) instead
       of iteratively with CGLS */
    bool directSolve = false;
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
// a total order, so that the line search does not depend on the sorting
// algorithm
inline bool operator<(const Delta& a, const Delta& b) {
  return a.delta < b.delta || (a.delta == b.delta && a.index < b.index);
}

/* Buffers of the L2_SVM_MFN solver. One workspace can serve several solves
   in turn (e.g. one per thread), which then do not allocate again unless
   they have more examples. */
struct SolverWorkspace {
  /* makes the buffers large enough for m examples with n values each; z and
     q are only needed by CGLS. Their contents are not kept. */
  void prepare(int m, int n, bool cgls);
  std::vector<double> z, q;  /* m, for CGLS */
  std::vector<double> r, p;  /* n */
  std::vector<Delta> deltas;  /* m, for the line search */
  std::vector<double> gram;  /* n x n, for the direct solve */
  std::vector<int> activeSubset;  /* m */
  std::vector<double> weightsBar;  /* n */
  std::vector<double> outputsBar;  /* m */
};

/* svmlin algorithms and their subroutines */

/* Conjugate Gradient for Sparse Linear Least Squares Problems */
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_{i in Subset} Data->C[i] (Y[i]- w' x_i)^2 */
/* over the active examples x_i, i in J[0..active) */
int CGLS(const AlgIn& set, const double lambda, const int cgitermax,
         const double epsilon, const int* J, int active,
         double* beta, double* o,
         double cpos, double cneg, SolverWorkspace& workspace);

/* Linear Modified Finite Newton L2-SVM*/
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_i Data->C[i] max(0,1 - Y[i] w' x_i)^2 */
/* Starts from Weights and their Outputs (the first set.m values are used),
   with the buffers of workspace */
int L2_SVM_MFN(const AlgIn& set, const options& Options,
               std::vector<double>& Weights,
               std::vector<double>& Outputs, double cpos, double cneg,
               SolverWorkspace& workspace);
/* the same with buffers of its own */
int L2_SVM_MFN(const AlgIn& set, const options& Options,
               std::vector<double>& Weights,
               std::vector<double>& Outputs, double cpos, double cneg);
double line_search(double* w, double* w_bar, double lambda, double* o,
                         double* o_bar, const double* Y, int d, int l,
                          double cpos, double cneg, Delta* deltas);
#endif
