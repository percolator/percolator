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

#ifndef SEQUENTIALSORT_H_
#define SEQUENTIALSORT_H_

#include <algorithm>
#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

/*
 * std::sort and std::nth_element on the calling thread, for calls inside
 * OpenMP parallel regions. With GCC's parallel mode (-D_GLIBCXX_PARALLEL)
 * the plain calls would start a nested team, which has one thread there,
 * for every call. Once the build no longer uses parallel mode these are the
 * plain calls.
 */
template <class RandomIt, class Compare>
inline void sequentialSort(RandomIt first, RandomIt last, Compare comp) {
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template <class RandomIt>
inline void sequentialSort(RandomIt first, RandomIt last) {
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

template <class RandomIt>
inline void sequentialNthElement(RandomIt first, RandomIt nth, RandomIt last) {
#ifdef _GLIBCXX_PARALLEL
  std::nth_element(first, nth, last, __gnu_parallel::sequential_tag());
#else
  std::nth_element(first, nth, last);
#endif
}

#endif /* SEQUENTIALSORT_H_ */
