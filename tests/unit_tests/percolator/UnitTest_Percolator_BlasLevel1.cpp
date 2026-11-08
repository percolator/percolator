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

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "BlasLevel1.h"

namespace {

/*
 * Reference: the increment-1 code of ddot_, daxpy_ and dscal_ of the former
 * src/dblas library (LINPACK, J. Dongarra), transcribed with its loop
 * unrolling and long indices. The solver's results depend on the summation
 * order of ddot.
 */
double referenceDdot(int n, const double* sx, const double* sy) {
  long nn = n;
  double stemp = 0.0;
  if (nn > 0) {
    long i;
    long m = nn - 4;
    for (i = 0; i < m; i += 5)
      stemp += sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * sy[i + 2] +
               sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];
    for (; i < nn; i++) stemp += sx[i] * sy[i];
  }
  return stemp;
}

void referenceDaxpy(int n, double ssa, const double* sx, double* sy) {
  long nn = n;
  if (nn > 0 && ssa != 0.0) {
    long i;
    long m = nn - 3;
    for (i = 0; i < m; i += 4) {
      sy[i] += ssa * sx[i];
      sy[i + 1] += ssa * sx[i + 1];
      sy[i + 2] += ssa * sx[i + 2];
      sy[i + 3] += ssa * sx[i + 3];
    }
    for (; i < nn; ++i) sy[i] += ssa * sx[i];
  }
}

void referenceDscal(int n, double ssa, double* sx) {
  long nn = n;
  if (nn > 0) {
    long i;
    long m = nn - 4;
    for (i = 0; i < m; i += 5) {
      sx[i] = ssa * sx[i];
      sx[i + 1] = ssa * sx[i + 1];
      sx[i + 2] = ssa * sx[i + 2];
      sx[i + 3] = ssa * sx[i + 3];
      sx[i + 4] = ssa * sx[i + 4];
    }
    for (; i < nn; ++i) sx[i] = ssa * sx[i];
  }
}

bool sameBits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() &&
         (a.empty() ||
          std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

// values over many orders of magnitude and of both signs, so that the
// summation order changes the rounding
std::vector<double> randomVector(std::mt19937& rng, int n) {
  std::uniform_real_distribution<double> mantissa(-1.0, 1.0);
  std::uniform_int_distribution<int> exponent(-20, 20);
  std::vector<double> v(static_cast<size_t>(n));
  for (double& x : v) x = std::ldexp(mantissa(rng), exponent(rng));
  return v;
}

}  // namespace

// Same bits as the reference for every length up to 70, which covers all
// remainders of the unrolled loops and the usual lengths of the solver's
// weight vectors (features plus one).
TEST(BlasLevel1Test, SameAsReference) {
  std::mt19937 rng(17);
  int numOrderSensitive = 0;
  for (int rep = 0; rep < 20; ++rep) {
    for (int n = 0; n <= 70; ++n) {
      std::vector<double> x = randomVector(rng, n);
      std::vector<double> y = randomVector(rng, n);
      const double expected = referenceDdot(n, x.data(), y.data());
      ASSERT_TRUE(sameBits(expected, blas::dot(n, x.data(), y.data())))
          << "n=" << n;
      // the plain left-to-right sum differs for some inputs, so the test
      // would notice a change of the order
      double plain = 0.0;
      for (int i = 0; i < n; ++i) plain += x[i] * y[i];
      if (!sameBits(plain, expected)) ++numOrderSensitive;

      for (double alpha : {0.75, -3.0e-7, 1.0, -1.0, 0.0}) {
        std::vector<double> yRef = y, yNew = y;
        referenceDaxpy(n, alpha, x.data(), yRef.data());
        blas::axpy(n, alpha, x.data(), yNew.data());
        ASSERT_TRUE(sameBits(yRef, yNew)) << "axpy n=" << n;
        std::vector<double> xRef = x, xNew = x;
        referenceDscal(n, alpha, xRef.data());
        blas::scal(n, alpha, xNew.data());
        ASSERT_TRUE(sameBits(xRef, xNew)) << "scal n=" << n;
      }
    }
  }
  EXPECT_GT(numOrderSensitive, 100);
}

// axpy with alpha 0 leaves y unchanged even if x holds a NaN or an infinity;
// scal by 0 does not skip them.
TEST(BlasLevel1Test, ZeroAlpha) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> x = {1.0, nan, inf, -2.0};
  std::vector<double> y = {3.0, 4.0, 5.0, 6.0};
  blas::axpy(4, 0.0, x.data(), y.data());
  EXPECT_TRUE(sameBits(std::vector<double>({3.0, 4.0, 5.0, 6.0}), y));
  blas::scal(4, 0.0, x.data());
  EXPECT_EQ(0.0, x[0]);
  EXPECT_TRUE(std::isnan(x[1]));
  EXPECT_TRUE(std::isnan(x[2]));
  EXPECT_EQ(0.0, x[3]);
  EXPECT_EQ(0.0, blas::dot(0, x.data(), y.data()));
}
