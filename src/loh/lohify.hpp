#ifndef _LOHIFY_HPP
#define _LOHIFY_HPP

#include <functional>
#include <algorithm>

#include "PrimitiveVector.hpp"

template <typename T>
T* triple_partition_naive(T*begin, T*end, T pivot) {
  auto geq_iterator = std::partition(begin, end, [pivot](auto val){return val>pivot;});
  std::partition(geq_iterator, end, [pivot](auto val){return val==pivot;});
  return geq_iterator;
}


// expected to behave roughly like alpha=2
template <typename T, int MINIMUM_LAYER_SIZE=128, bool RANDOMIZE=true, bool FORCE_FIRST_LAYER_TO_HAVE_SIZE_1=true>
void max_quick_lohify(T*__restrict x, T*__restrict x_end, PrimitiveVector<unsigned long> & partition_ranks, std::function<bool(const T&,const T&)> compare) {
  unsigned long n = x_end - x;
  while (n > MINIMUM_LAYER_SIZE) {
    // base case check
    unsigned long pivot_index;
    if (RANDOMIZE)
      pivot_index = rand() % n;
    else
      pivot_index = n/2;

    T pivot_val = x[pivot_index];
    std::swap(x[n-1], x[pivot_index]);

    T*second_half_of_x = triple_partition_naive<T>(x, x_end, pivot_val);
    
    // pivot element will be moved to last element of the right half; swap it back so that it begins the right half:
    std::swap(*second_half_of_x, *(x_end-1));

    if (second_half_of_x > x)
      partition_ranks.push_back(second_half_of_x - x);

    x_end = second_half_of_x;
    n = second_half_of_x - x;
  }

  if (FORCE_FIRST_LAYER_TO_HAVE_SIZE_1) {
    std::sort(x,x_end,[](auto lhs, auto rhs) { return lhs > rhs; });

    for (auto iter=x_end-1; iter>x; --iter)
      partition_ranks.push_back(iter-x);
  }

  std::reverse(partition_ranks.begin(), partition_ranks.end());
}

#endif
