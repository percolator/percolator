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
 * Unit tests for parallelSort: for the total orders of the PSM sorts it must
 * give exactly the order of std::sort, on any number of threads.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ParallelSort.h"
#include "PsmStore.h"
#include "ScoreHolder.h"

namespace {

// Runs body with each of the thread counts, and restores the thread count.
template <class F>
void forThreadCounts(const F& body) {
#ifdef _OPENMP
  const int saved = omp_get_max_threads();
  for (int threads : {1, 2, 3, 4, 5, 8}) {
    omp_set_num_threads(threads);
    body(threads);
  }
  omp_set_num_threads(saved);
#else
  body(1);
#endif
}

unsigned int nextRandom(unsigned int& seed) {
  seed = seed * 1103515245u + 12345u;
  return seed >> 8;
}

// PSMs with coarse scores, scans, masses, labels and ids, so that the
// comparators often reach their last tie-breaks; the ids repeat, as with
// duplicated input rows, and only the line numbers tell those apart.
class PsmSet {
 public:
  explicit PsmSet(int n) : ids_(static_cast<size_t>(n)),
      peptides_(static_cast<size_t>(n)), psms_(static_cast<size_t>(n)) {
    unsigned int seed = 4711u;
    const char* aminoAcids = "ACDE";
    for (int i = 0; i < n; ++i) {
      ids_[i] = "psm" + std::to_string(nextRandom(seed) % (n / 4 + 1));
      std::string peptide = "K.";
      for (unsigned int j = 0, len = 1 + nextRandom(seed) % 3; j < len; ++j) {
        peptide += aminoAcids[nextRandom(seed) % 4];
      }
      peptides_[i] = peptide + ".R";
    }
    for (int i = 0; i < n; ++i) {
      PSMDescription& psm = psms_[i];
      store_.storeStrings(psm, ids_[i], peptides_[i], {});
      psm.setLineNr(static_cast<unsigned int>(i));
      psm.scan = nextRandom(seed) % (n / 3 + 1);
      psm.specFileNr = nextRandom(seed) % 2;
      psm.expMass = 1000.0 + nextRandom(seed) % 3;
      const LabelType label =
          nextRandom(seed) % 2 ? LabelType::TARGET : LabelType::DECOY;
      scores_.push_back(ScoreHolder(
          static_cast<double>(nextRandom(seed) % 50) / 4.0, label, &psm));
    }
  }
  const std::vector<ScoreHolder>& scores() const { return scores_; }

 private:
  std::vector<std::string> ids_, peptides_;
  PsmStore store_;  // holds the ids and peptides of psms_
  std::vector<PSMDescription> psms_;
  std::vector<ScoreHolder> scores_;
};

template <class Compare>
void expectSameOrderAsStdSort(const std::vector<ScoreHolder>& input,
                              Compare comp, const char* name) {
  std::vector<ScoreHolder> expected = input;
  std::sort(expected.begin(), expected.end(), comp);
  forThreadCounts([&](int threads) {
    std::vector<ScoreHolder> sorted = input;
    parallelSort(sorted.begin(), sorted.end(), comp);
    ASSERT_EQ(expected.size(), sorted.size());
    for (size_t i = 0; i < sorted.size(); ++i) {
      ASSERT_EQ(expected[i].pPSM, sorted[i].pPSM)
          << name << " with " << threads << " threads, position " << i;
    }
  });
}

}  // namespace

// The comparators of the PSM sorts on random PSMs with many ties, which
// only the PSM id and the line number break.
TEST(ParallelSortTest, ScoreHolderOrdersEqualStdSort) {
  PsmSet set(120000);
  std::vector<ScoreHolder> input = set.scores();
  expectSameOrderAsStdSort(input, std::greater<ScoreHolder>(), "greater");
  expectSameOrderAsStdSort(input, OrderScanHash(), "OrderScanHash");
  expectSameOrderAsStdSort(input, lexicOrderProb(), "lexicOrderProb");
  expectSameOrderAsStdSort(input, OrderScanMassCharge(), "OrderScanMassCharge");
  // already sorted, reversed, and a concatenation of sorted bins, as in
  // Scores::merge
  std::sort(input.begin(), input.end(), std::greater<ScoreHolder>());
  expectSameOrderAsStdSort(input, std::greater<ScoreHolder>(), "sorted");
  std::reverse(input.begin(), input.end());
  expectSameOrderAsStdSort(input, std::greater<ScoreHolder>(), "reversed");
  input = set.scores();
  for (size_t bin = 0; bin < 3; ++bin) {
    std::sort(input.begin() + input.size() * bin / 3,
              input.begin() + input.size() * (bin + 1) / 3,
              std::greater<ScoreHolder>());
  }
  expectSameOrderAsStdSort(input, std::greater<ScoreHolder>(), "bins");
  // all scores tied, so that the split elements are chosen among ties
  for (ScoreHolder& sh : input) sh.score = 1.0;
  expectSameOrderAsStdSort(input, std::greater<ScoreHolder>(), "tied");
}

// Values with many duplicates and lengths around the parallel threshold.
TEST(ParallelSortTest, ValuesEqualStdSort) {
  for (int n : {0, 1, 2, 1000, 32767, 32768, 32769, 65536 + 17, 300001}) {
    unsigned int seed = static_cast<unsigned int>(n);
    std::vector<std::pair<double, bool> > input;
    for (int i = 0; i < n; ++i) {
      const unsigned int r = nextRandom(seed);
      input.push_back(std::make_pair(static_cast<double>(r % 1000), r % 3 == 0));
    }
    std::vector<std::pair<double, bool> > expected = input;
    std::sort(expected.begin(), expected.end(),
              std::greater<std::pair<double, bool> >());
    forThreadCounts([&](int threads) {
      std::vector<std::pair<double, bool> > sorted = input;
      parallelSort(sorted.begin(), sorted.end(),
                   std::greater<std::pair<double, bool> >());
      EXPECT_TRUE(sorted == expected) << n << " values, " << threads
                                      << " threads";
    });
  }
}

// Inside a parallel region it sorts on the calling thread.
TEST(ParallelSortTest, SortsInsideParallelRegion) {
  std::vector<std::vector<int> > inputs(4);
  unsigned int seed = 99u;
  for (std::vector<int>& v : inputs) {
    for (int i = 0; i < 100000; ++i) {
      v.push_back(static_cast<int>(nextRandom(seed) % 5000));
    }
  }
  std::vector<std::vector<int> > sorted = inputs;
#pragma omp parallel for num_threads(2)
  for (int i = 0; i < 4; ++i) {
    parallelSort(sorted[i].begin(), sorted[i].end(), std::less<int>());
  }
  for (int i = 0; i < 4; ++i) {
    std::sort(inputs[i].begin(), inputs[i].end());
    EXPECT_TRUE(inputs[i] == sorted[i]) << "vector " << i;
  }
}

// An exception thrown by the comparator on another thread than the calling
// one reaches the caller instead of terminating the program.
TEST(ParallelSortTest, PassesOnComparatorExceptions) {
  std::vector<int> input(200000);
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = static_cast<int>((i * 7919u) % input.size());
  }
  auto failsOnSecondThread = [](int a, int b) {
#ifdef _OPENMP
    if (omp_get_thread_num() == 1) throw std::runtime_error("comparator");
#endif
    return a < b;
  };
  forThreadCounts([&](int threads) {
    std::vector<int> values = input;
    if (threads > 1) {
      EXPECT_THROW(
          parallelSort(values.begin(), values.end(), failsOnSecondThread),
          std::runtime_error) << threads << " threads";
    } else {
      parallelSort(values.begin(), values.end(), failsOnSecondThread);
      EXPECT_TRUE(std::is_sorted(values.begin(), values.end()));
    }
  });
}
