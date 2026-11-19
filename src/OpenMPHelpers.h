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

#endif  // OPENMP_HELPERS_H_
