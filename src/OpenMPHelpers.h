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

#ifndef OPENMP_HELPERS_H_
#define OPENMP_HELPERS_H_

#include <atomic>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

// Whether the calling thread runs in an active parallel region; false
// without OpenMP. For the clauses of OpenMP directives: Clang parses them
// also with -fopenmp-simd, which does not define _OPENMP, so that omp.h is
// not included and omp_in_parallel() is not declared there.
inline bool inParallel() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

// An exception that leaves an OpenMP parallel region ends the process
// (std::terminate), also in a program that calls Percolator as a library
// and would catch it. Each iteration of a parallel loop therefore catches
// its own exception and passes it to capture(); the first one is kept,
// the iterations that start later do nothing, and rethrow() throws the
// kept exception on the calling thread after the region:
//
//   ParallelErrors errors;
//   #pragma omp parallel for
//   for (int i = 0; i < n; ++i) {
//     if (errors.failed()) continue;
//     try {
//       ...
//     } catch (...) {
//       errors.capture();
//     }
//   }
//   errors.rethrow();
//
// Only a critical section is used, so that it works with OpenMP 2.0 (MSVC).
class ParallelErrors {
 public:
  ParallelErrors() : failed_(false) {}
  ParallelErrors(const ParallelErrors&) = delete;
  ParallelErrors& operator=(const ParallelErrors&) = delete;

  // whether an iteration has failed; the others may then skip their work
  bool failed() const { return failed_.load(std::memory_order_relaxed); }

  // to be called in a catch block: keeps the exception being handled if it
  // is the first
  void capture() {
#pragma omp critical(percolator_parallel_errors)
    {
      if (!error_) {
        error_ = std::current_exception();
      }
    }
    failed_.store(true, std::memory_order_relaxed);
  }

  // after the parallel region: throws the kept exception, if any
  void rethrow() const {
    if (error_) {
      std::rethrow_exception(error_);
    }
  }

 private:
  std::exception_ptr error_;
  std::atomic<bool> failed_;
};

#endif  // OPENMP_HELPERS_H_
