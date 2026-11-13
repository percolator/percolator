/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
/*
 * The three level-1 BLAS operations of the SVM solver (ssl.cpp), for
 * contiguous vectors. They replace ddot_, daxpy_ and dscal_ of the former
 * src/dblas library and keep its loops (LINPACK, J. Dongarra), so that they
 * do the same floating-point operations in the same order, and the solver's
 * results stay the same to the last bit, with as many instructions.
 */
#ifndef BLAS_LEVEL1_H_
#define BLAS_LEVEL1_H_

namespace blas {

/**
 * x'y. The products are added in groups of five, each group summed from left
 * to right and then added to the running sum, and the remaining products one
 * by one: the order of the reference ddot, on which the result depends.
 */
inline double dot(int n, const double* x, const double* y) {
  const long nn = n;
  double sum = 0.0;
  long i = 0;
  for (; i < nn - 4; i += 5) {
    sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] +
           x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
  }
  for (; i < nn; ++i) {
    sum += x[i] * y[i];
  }
  return sum;
}

/** y += alpha * x. Like daxpy, does nothing if alpha is 0, even where x holds
    an infinity or a NaN. */
inline void axpy(int n, double alpha, const double* x, double* y) {
  if (alpha == 0.0) return;
  const long nn = n;
  long i = 0;
  for (; i < nn - 3; i += 4) {
    y[i] += alpha * x[i];
    y[i + 1] += alpha * x[i + 1];
    y[i + 2] += alpha * x[i + 2];
    y[i + 3] += alpha * x[i + 3];
  }
  for (; i < nn; ++i) {
    y[i] += alpha * x[i];
  }
}

/** x = alpha * x */
inline void scal(int n, double alpha, double* x) {
  const long nn = n;
  long i = 0;
  for (; i < nn - 4; i += 5) {
    x[i] = alpha * x[i];
    x[i + 1] = alpha * x[i + 1];
    x[i + 2] = alpha * x[i + 2];
    x[i + 3] = alpha * x[i + 3];
    x[i + 4] = alpha * x[i + 4];
  }
  for (; i < nn; ++i) {
    x[i] = alpha * x[i];
  }
}

}  // namespace blas

#endif /* BLAS_LEVEL1_H_ */
