#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <cmath>
#include <algorithm>

namespace phos_loc {

template<typename Iterator>
bool NextCombination(const Iterator first, Iterator k, const Iterator last);
unsigned int NChooseK(unsigned int n, unsigned int k);
double LogNChooseK(unsigned int n, unsigned int k);
double BinomialPMF(unsigned int num_trials,
                    unsigned int num_successes,
                    double prob_success);
double CumulativeBinomialProbability(unsigned int num_trials,
                                      unsigned int num_successes,
                                      double prob_success);
} // namespace phos_loc

#endif // UTILS_H_
