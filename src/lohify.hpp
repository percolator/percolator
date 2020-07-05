#ifndef _LOHIFY_HPP
#define _LOHIFY_HPP

#include <functional>
#include <algorithm>

#include "PrimitiveVector.hpp"

// todo:  SIZE_TO_SORT should be more clever, looking at density of pivots needed rather than total elements
template <typename T, int SIZE_TO_SORT=16>
void max_random_lohify(T*x, T*x_end, unsigned long*partition_ranks, unsigned long*partition_ranks_end, unsigned long n_eliminated_from_left=0) {
    if (x >= x_end || partition_ranks >= partition_ranks_end)
      return;

    unsigned long n=x_end-x;
    if (n <= SIZE_TO_SORT) {
      std::sort(x,x_end,[](auto lhs, auto rhs){return lhs > rhs;});
      return;
    }

    // Note: more careful choice of pivot_index could result in optimal runtime (as per median-of-medians)
    unsigned long pivot_index = rand() % n;
    T pivot_val = x[pivot_index];

    std::swap(x[n-1], x[pivot_index]);
    T*second_half_of_x = std::partition(x, x_end, [pivot_val](auto val) {return val > pivot_val;});
    // pivot element will be moved to last element of the right half; swap it back so that it begins the right half:
    std::swap(*second_half_of_x, *(x_end-1));

    unsigned long n_left = second_half_of_x - x;
    // Do not bother partitioning partition_ranks, because it is already
    // sorted; instead, search in log2(n_partition_ranks) time:
    unsigned long partition_index_to_split = n_left + n_eliminated_from_left;
    unsigned long* second_half_of_partition_ranks = std::lower_bound(partition_ranks, partition_ranks_end, partition_index_to_split);

    if (x < second_half_of_x && partition_ranks < second_half_of_partition_ranks)
      max_random_lohify(x, second_half_of_x, partition_ranks, second_half_of_partition_ranks, n_eliminated_from_left);

    // Add second_half_of_x+1 because pivot is already in sorted position. So it should
    // be considered by neither left nor right recursion:

    if (second_half_of_partition_ranks < partition_ranks_end && n_left+n_eliminated_from_left == *second_half_of_partition_ranks) {
      // if lower bound is one of the partitions we want (you fish your wish!), remove it from the list of partitions.
      if (second_half_of_x+1 < x_end && second_half_of_partition_ranks+1 < partition_ranks_end)
	max_random_lohify(second_half_of_x+1, x_end, second_half_of_partition_ranks+1, partition_ranks_end, n_eliminated_from_left+n_left+1);
    }  
    else
      if (second_half_of_x+1 < x_end && second_half_of_partition_ranks < partition_ranks_end)
	max_random_lohify(second_half_of_x+1, x_end, second_half_of_partition_ranks, partition_ranks_end, n_eliminated_from_left+n_left+1);

}

// expected to behave roughly like alpha=2
template <typename T, int MINIMUM_LAYER_SIZE=128, bool RANDOMIZE=true, bool FORCE_FIRST_LAYER_TO_HAVE_SIZE_1=true>
void max_quick_lohify(T*__restrict x, T*__restrict x_end, PrimitiveVector<unsigned long> & partition_ranks, std::function<bool(const T&,const T&)> compare) {
  unsigned long n = x_end - x;
  while (n > MINIMUM_LAYER_SIZE) {

    unsigned long pivot_index;
    if (RANDOMIZE)
      pivot_index = rand() % n;
    else
      pivot_index = n/2;
    T pivot_val = x[pivot_index];

    std::swap(x[n-1], x[pivot_index]);

    T*second_half_of_x = std::partition(x, x_end, [pivot_val,&compare](auto val) {return compare(val,pivot_val);});
    // pivot element will be moved to last element of the right half; swap it back so that it begins the right half:
    std::swap(*second_half_of_x, *(x_end-1));
    partition_ranks.push_back(second_half_of_x - x);

    x_end = second_half_of_x;
    n = second_half_of_x - x;
  }

  if (FORCE_FIRST_LAYER_TO_HAVE_SIZE_1) {
    std::sort(x,x_end,compare);

    for (auto iter=x_end-1; iter>x; --iter)
      partition_ranks.push_back(iter-x);
  }

  std::reverse(partition_ranks.begin(), partition_ranks.end());
}




#endif
