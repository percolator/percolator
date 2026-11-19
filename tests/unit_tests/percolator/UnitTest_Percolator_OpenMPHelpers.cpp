/*******************************************************************************
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
 * Unit tests for OpenMPHelpers.h: ParallelErrors passes an exception of an
 * iteration of a parallel loop on to the calling thread, where it can be
 * caught; without it the process ends (std::terminate).
 */

#include <gtest/gtest.h>

#include <new>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "OpenMPHelpers.h"

namespace {

int threadNumber() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

int teamSize() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

}  // namespace

// The exception of one iteration reaches the calling thread after the loop,
// whichever iteration and thread it came from, and the iteration that threw
// did not finish.
TEST(ParallelErrorsTest, RethrowsTheExceptionOfOneIteration) {
  const int n = 1000;
  for (int failing : {0, 1, 499, n - 1}) {
    ParallelErrors errors;
    std::vector<int> finished(n, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
    for (int i = 0; i < n; ++i) {
      if (errors.failed()) {
        continue;
      }
      try {
        if (i == failing) {
          throw std::runtime_error("iteration " + std::to_string(i));
        }
        finished[i] = 1;
      } catch (...) {
        errors.capture();
      }
    }
    EXPECT_TRUE(errors.failed());
    bool caught = false;
    try {
      errors.rethrow();
    } catch (const std::runtime_error& e) {
      caught = true;
      EXPECT_EQ("iteration " + std::to_string(failing), e.what());
    }
    EXPECT_TRUE(caught) << "iteration " << failing;
    EXPECT_EQ(0, finished[failing]);
  }
}

// On one thread the iterations run in order: those before the failing one
// finish, and those after it do nothing.
TEST(ParallelErrorsTest, SkipsTheIterationsAfterTheError) {
  const int n = 100;
  const int failing = 40;
  ParallelErrors errors;
  std::vector<int> finished(n, 0);
#pragma omp parallel for schedule(static) num_threads(1)
  for (int i = 0; i < n; ++i) {
    if (errors.failed()) {
      continue;
    }
    try {
      if (i == failing) {
        throw std::runtime_error("iteration " + std::to_string(i));
      }
      finished[i] = 1;
    } catch (...) {
      errors.capture();
    }
  }
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(i < failing ? 1 : 0, finished[i]) << "iteration " << i;
  }
  EXPECT_THROW(errors.rethrow(), std::runtime_error);
}

// An exception thrown on another thread than the calling one keeps its
// type.
TEST(ParallelErrorsTest, RethrowsAnExceptionOfAnotherThread) {
  const int n = 64;
  ParallelErrors errors;
  int threads = 1;
#pragma omp parallel num_threads(4)
  {
#pragma omp single
    threads = teamSize();
#pragma omp for schedule(static, 1)
    for (int i = 0; i < n; ++i) {
      if (errors.failed()) {
        continue;
      }
      try {
        if (threadNumber() == 1) {
          throw std::bad_alloc();
        }
      } catch (...) {
        errors.capture();
      }
    }
  }
  if (threads > 1) {
    EXPECT_TRUE(errors.failed());
    EXPECT_THROW(errors.rethrow(), std::bad_alloc);
  } else {
    EXPECT_FALSE(errors.failed());
    EXPECT_NO_THROW(errors.rethrow());
  }
}

// Of several exceptions, the first one captured is rethrown. The loop does
// not check failed(), so that every iteration throws, and the order of the
// captures is recorded in a critical section of its own around capture().
TEST(ParallelErrorsTest, KeepsTheFirstException) {
  const int n = 64;
  ParallelErrors errors;
  std::vector<int> order;
#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
  for (int i = 0; i < n; ++i) {
    try {
      throw std::runtime_error("iteration " + std::to_string(i));
    } catch (...) {
#pragma omp critical(parallel_errors_test_order)
      {
        order.push_back(i);
        errors.capture();
      }
    }
  }
  ASSERT_EQ(static_cast<size_t>(n), order.size());
  bool caught = false;
  try {
    errors.rethrow();
  } catch (const std::runtime_error& e) {
    caught = true;
    EXPECT_EQ("iteration " + std::to_string(order.front()), e.what());
  }
  EXPECT_TRUE(caught);
}

// Without an exception every iteration runs and rethrow does nothing.
TEST(ParallelErrorsTest, NothingToRethrowWithoutAnException) {
  const int n = 1000;
  ParallelErrors errors;
  std::vector<int> finished(n, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
  for (int i = 0; i < n; ++i) {
    if (errors.failed()) {
      continue;
    }
    try {
      finished[i] = 1;
    } catch (...) {
      errors.capture();
    }
  }
  EXPECT_FALSE(errors.failed());
  EXPECT_NO_THROW(errors.rethrow());
  EXPECT_EQ(std::vector<int>(n, 1), finished);
}

// inParallel() is true in a parallel region of more than one thread only.
TEST(InParallelTest, OnlyInsideAnActiveParallelRegion) {
  EXPECT_FALSE(inParallel());
  int threads = 1;
  bool inside = false;
#pragma omp parallel num_threads(2)
  {
#pragma omp single
    {
      threads = teamSize();
      inside = inParallel();
    }
  }
  EXPECT_EQ(threads > 1, inside);
  EXPECT_FALSE(inParallel());
}
