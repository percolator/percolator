#include "utils.h"

namespace utils {


template bool NextCombination<std::vector<int>::iterator>(
    const std::vector<int>::iterator first,
    std::vector<int>::iterator k,
    const std::vector<int>::iterator last);

template<typename Iterator>
bool NextCombination(const Iterator first, Iterator k, const Iterator last) {
  /* Credits: Thomas Draper */
  if ((first == last) || (first == k) || (last == k))
    return false;

  Iterator itr1 = first;
  Iterator itr2 = last;
  ++itr1;
  if (last == itr1)
    return false;

  itr1 = last;
  --itr1;
  itr1 = k;
  --itr2;
  while (first != itr1) {
    if (*--itr1 < *itr2) {
      Iterator j = k;
      while (!(*itr1 < *j)) ++j;
      std::iter_swap(itr1,j);
      ++itr1;
      ++j;
      itr2 = k;
      std::rotate(itr1,j,last);
      while (last != j) {
        ++j;
        ++itr2;
      }
      std::rotate(k,itr2,last);
      return true;
    }
  }
  std::rotate(first,k,last);
  return false;
}

unsigned int NChooseK(unsigned int n, unsigned int k) {
  if (k > n) return 0;
  if (k * 2 > n) k = n - k;
  if (k == 0) return 1;

  unsigned int result = n;
  for (unsigned int i = 2; i <= k; ++i) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

double BinomialPMF(unsigned int num_trials,
                   unsigned int num_successes,
                   double prob_success) {
  return NChooseK(num_trials, num_successes) * pow(prob_success, num_successes) *
      pow(1 - prob_success, num_trials - num_successes);
}

// Pr(X>=x)
double CumulativeBinomialProbability(unsigned int num_trials,
                                     unsigned int num_successes,
                                     double prob_success) {
  if (num_successes == 0)
    return 1.0;
  double cum_prob = 0;
  for (unsigned int k = num_successes; k <= num_trials; ++k) {
    cum_prob += BinomialPMF(num_trials, k, prob_success);
  }
  return cum_prob;
}

} // namespace utils
