/*
 * Unit tests for the memory reporting of the stage lines and for the
 * memory handling of the SVM training.
 */

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "PSMDescription.h"
#include "Scores.h"
#include "Timer.h"
#include "ssl.h"

TEST(MemoryUsageTest, parsesProcStatus) {
  std::istringstream status(
      "Name:\tpercolator\n"
      "VmPeak:\t 2000000 kB\n"
      "VmHWM:\t 1050620 kB\n"
      "VmRSS:\t  998400 kB\n"
      "RssAnon:\t  990000 kB\n"
      "Threads:\t3\n");
  MemoryUsage usage;
  parseProcStatus(status, usage);
  EXPECT_EQ(998400, usage.rssKb);
  EXPECT_EQ(1050620, usage.peakRssKb);
  EXPECT_EQ(-1, usage.heapKb);
}

TEST(MemoryUsageTest, keepsMissingAndMalformedFieldsUnset) {
  std::istringstream status("Name:\tx\nVmRSS:\tunknown\n");
  MemoryUsage usage;
  parseProcStatus(status, usage);
  EXPECT_EQ(-1, usage.rssKb);
  EXPECT_EQ(-1, usage.peakRssKb);
}

TEST(MemoryUsageTest, formatsAvailableFieldsInMiB) {
  MemoryUsage usage;
  EXPECT_EQ("", formatMemoryUsage(usage));
  usage.rssKb = 998400;
  usage.peakRssKb = 1050620;
  EXPECT_EQ(", rss 975 MiB, peak 1025 MiB", formatMemoryUsage(usage));
  usage.heapKb = 2047;
  EXPECT_EQ(", rss 975 MiB, peak 1025 MiB, heap 1 MiB",
            formatMemoryUsage(usage));
}

#if defined(__linux__)
TEST(MemoryUsageTest, readsOwnProcess) {
  MemoryUsage before = getMemoryUsage();
  ASSERT_GT(before.rssKb, 0);
  EXPECT_GE(before.peakRssKb, before.rssKb);
  // touching 64 MiB raises the resident set size and its peak
  std::vector<char> block(64u << 20, 1);
  MemoryUsage after = getMemoryUsage();
  EXPECT_GE(after.rssKb, before.rssKb + (48 << 10));
  EXPECT_GE(after.peakRssKb, after.rssKb);
  EXPECT_EQ(1, block[block.size() / 2]);
}
#endif

// merge releases each cross-validation bin after copying it, and gives the
// same merged list as merging copies of the bins
TEST(TrainingMemoryTest, mergeEmptiesTheBins) {
  FeatureNames::resetNumFeatures();
  FeatureNames::setNumFeatures(1);
  const int numPsms = 600;
  std::vector<PSMDescription> psms(numPsms);
  std::vector<Scores> bins(3, Scores(false));
  unsigned int seed = 4711u;
  for (int i = 0; i < numPsms; ++i) {
    seed = seed * 1103515245u + 12345u;
    const bool isTarget = (seed >> 16) % 2 == 0;
    psms[i].scan = static_cast<unsigned int>(i);
    psms[i].setLineNr(static_cast<unsigned int>(i));
    const double score =
        static_cast<double>((seed >> 4) % 50) / 10.0 + (isTarget ? 1.0 : 0.0);
    bins[i % 3].addScoreHolder(ScoreHolder(
        score, isTarget ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
  }
  for (Scores& bin : bins) {
    bin.recalculateSizes();
  }
  std::vector<Scores> copies = bins;
  std::vector<std::vector<double> > weights(3, std::vector<double>(2, 1.0));
  Scores merged(false), mergedCopies(false);
  merged.merge(bins, 0.05, true, weights);
  for (const Scores& bin : bins) {
    EXPECT_EQ(0u, bin.size());
    EXPECT_TRUE(bin.begin() == bin.end());
  }
  // the reference: the bins of copies are merged by hand
  std::vector<ScoreHolder> expected;
  for (Scores& copy : copies) {
    std::sort(copy.begin(), copy.end(), std::greater<ScoreHolder>());
    copy.calcQvals(0.05);
    expected.insert(expected.end(), copy.begin(), copy.end());
  }
  std::sort(expected.begin(), expected.end(), std::greater<ScoreHolder>());
  ASSERT_EQ(static_cast<unsigned int>(numPsms), merged.size());
  std::size_t ix = 0;
  for (const ScoreHolder& sh : merged) {
    EXPECT_EQ(expected[ix].pPSM, sh.pPSM);
    EXPECT_EQ(expected[ix].score, sh.score);
    EXPECT_EQ(expected[ix].q, sh.q);
    ++ix;
  }
  FeatureNames::resetNumFeatures();
}

// A training set of m examples with numFeatures features, two of them
// informative, and labels +1 or -1. The rows are kept in storage, with the
// intercept 1.0 after the features of each row.
static void makeTrainingSet(int m, int numFeatures, unsigned int seed,
                            std::vector<double>& storage, AlgIn& data) {
  const int n = numFeatures + 1;
  storage.assign(static_cast<std::size_t>(m) * n, 1.0);
  std::vector<double> labels(static_cast<std::size_t>(m));
  for (int i = 0; i < m; ++i) {
    seed = seed * 1103515245u + 12345u;
    const bool positive = (seed >> 16) % 3 != 0;
    labels[i] = positive ? 1.0 : -1.0;
    for (int k = 0; k < numFeatures; ++k) {
      seed = seed * 1103515245u + 12345u;
      double value = static_cast<double>((seed >> 8) % 1000) / 500.0 - 1.0;
      if (k < 2 && positive) value += 0.8;
      storage[static_cast<std::size_t>(i) * n + k] = value;
    }
  }
  data.clear();
  for (int i = 0; i < m; ++i) {  // negatives first, as Scores adds them
    if (labels[i] < 0) data.add(&storage[static_cast<std::size_t>(i) * n], -1);
  }
  data.negatives = data.m;
  for (int i = 0; i < m; ++i) {
    if (labels[i] > 0) data.add(&storage[static_cast<std::size_t>(i) * n], 1);
  }
  data.positives = data.m - data.negatives;
}

// Solves that share one workspace, with more and then fewer examples than
// the solve before, give the same weights and outputs as solves with their
// own buffers, with the direct solve and with CGLS.
TEST(TrainingMemoryTest, reusedSolverWorkspaceGivesTheSameSolution) {
  const int numFeatures = 6;
  const int n = numFeatures + 1;
  for (bool directSolve : {true, false}) {
    options opts;
    opts.lambda = 1.0;
    opts.lambda_u = 1.0;
    opts.epsilon = EPSILON;
    opts.cgitermax = CGITERMAX;
    opts.mfnitermax = MFNITERMAX;
    opts.directSolve = directSolve;
    SolverWorkspace shared;
    unsigned int seed = 17u;
    for (int m : {300, 700, 200, 701}) {
      std::vector<double> storage;
      AlgIn data(n);
      makeTrainingSet(m, numFeatures, seed++, storage, data);
      for (double cneg : {1.0, 3.0}) {
        std::vector<double> w1(n, 0.0), o1(m, 0.0);
        std::vector<double> w2(n, 0.0), o2(m + 5, 0.0);  // longer is allowed
        L2_SVM_MFN(data, opts, w1, o1, 1.0, cneg);
        L2_SVM_MFN(data, opts, w2, o2, 1.0, cneg, shared);
        for (int k = 0; k < n; ++k) {
          EXPECT_EQ(w1[k], w2[k]) << "m " << m << " k " << k;
        }
        for (int i = 0; i < m; ++i) {
          EXPECT_EQ(o1[i], o2[i]) << "m " << m << " i " << i;
        }
        EXPECT_NE(0.0, w1[0]);
      }
    }
    // the CGLS buffers are only allocated for CGLS
    EXPECT_EQ(directSolve ? 0u : 701u, shared.z.size());
    EXPECT_EQ(directSolve ? 0u : 701u, shared.q.size());
    EXPECT_EQ(701u, shared.deltas.size());
  }
}

// Rows of a pool set up for in-place training end with the intercept 1.0
// of the SVM, also when a row is used again after deallocate, and do not
// overlap; by default the rows hold the features only.
TEST(TrainingMemoryTest, featurePoolRowsForInPlaceTraining) {
  const unsigned int numFeatures = 30;
  for (bool inPlace : {false, true}) {
    FeatureMemoryPool pool;
    pool.setSvmRowsInPlace(inPlace);
    pool.createPool(numFeatures);
    EXPECT_EQ(inPlace, pool.hasInterceptColumn());
    EXPECT_EQ(inPlace, pool.canTrainInPlace(numFeatures));
    EXPECT_FALSE(pool.canTrainInPlace(numFeatures - 1));
    EXPECT_EQ(numFeatures + (inPlace ? 1u : 0u), pool.rowStride());
    // more rows than two blocks hold
    const unsigned int numRows = 2 * ((1u << 20) / pool.rowStride()) + 7;
    std::vector<double*> rows;
    for (unsigned int i = 0; i < numRows; ++i) {
      double* row = pool.allocate();
      if (inPlace) {
        EXPECT_EQ(1.0, row[numFeatures]);
      }
      for (unsigned int k = 0; k < numFeatures; ++k) {
        row[k] = static_cast<double>(i) + k / 100.0;
      }
      rows.push_back(row);
    }
    for (unsigned int i = 0; i < numRows; ++i) {
      ASSERT_EQ(static_cast<double>(i), rows[i][0]);
      ASSERT_EQ(static_cast<double>(i) + 0.29, rows[i][numFeatures - 1]);
      if (inPlace) {
        ASSERT_EQ(1.0, rows[i][numFeatures]);
      }
      if (i % ((1u << 20) / pool.rowStride()) != 0) {  // within a block
        ASSERT_EQ(rows[i - 1] + pool.rowStride(), rows[i]);
      }
    }
    if (inPlace) {
      // a used row gets its intercept back when it is allocated again
      rows[5][numFeatures] = -3.0;
      pool.deallocate(rows[5]);
      double* again = pool.allocate();
      EXPECT_EQ(rows[5], again);
      EXPECT_EQ(1.0, again[numFeatures]);
    }
    pool.destroyPool();
  }
}

// The training sets made by Scores hold the features of the PSMs and the
// intercept 1.0, as copies by default and as the rows themselves in place.
TEST(TrainingMemoryTest, trainingSetCopiesOrPointsToTheFeatureRows) {
  const unsigned int numFeatures = 3;
  FeatureNames::resetNumFeatures();
  FeatureNames::setNumFeatures(numFeatures);
  for (bool inPlace : {false, true}) {
    FeatureMemoryPool pool;
    pool.setSvmRowsInPlace(inPlace);
    pool.createPool(numFeatures);
    std::vector<PSMDescription> psms(8);
    Scores scores(false);
    for (unsigned int i = 0; i < psms.size(); ++i) {
      psms[i].features = pool.allocate();
      for (unsigned int k = 0; k < numFeatures; ++k) {
        psms[i].features[k] = i + 0.5 * k;
      }
      psms[i].scan = i;
      ScoreHolder sh(10.0 - i, i % 2 ? LabelType::DECOY : LabelType::TARGET,
                     &psms[i]);
      sh.q = i < 4 ? 0.0 : 0.5;
      scores.addScoreHolder(sh);
    }
    AlgIn data(static_cast<int>(numFeatures) + 1, inPlace);
    scores.generateNegativeTrainingSet(data, 1.0);
    scores.generatePositiveTrainingSet(data, 0.01, 1.0, false);
    ASSERT_EQ(4, data.negatives);  // the decoys 1, 3, 5, 7
    ASSERT_EQ(2, data.positives);  // the targets 0, 2 with q <= 0.01
    ASSERT_EQ(6, data.m);
    const unsigned int expected[] = {1, 3, 5, 7, 0, 2};
    for (int i = 0; i < data.m; ++i) {
      const double* features = psms[expected[i]].features;
      if (inPlace) {
        EXPECT_EQ(features, data.row(i));
      } else {
        EXPECT_NE(features, data.row(i));
        if (i > 0) {
          EXPECT_EQ(data.row(i - 1) + data.n, data.row(i));  // contiguous
        }
      }
      for (unsigned int k = 0; k < numFeatures; ++k) {
        EXPECT_EQ(features[k], data.row(i)[k]);
      }
      EXPECT_EQ(1.0, data.row(i)[numFeatures]);
      EXPECT_EQ(i < 4 ? -1.0 : 1.0, data.Y[i]);
    }
    pool.destroyPool();
  }
  FeatureNames::resetNumFeatures();
}

// Copied rows stay valid when the matrix grows without a reserve and when
// the training set is copied or moved.
TEST(TrainingMemoryTest, copiedRowsFollowTheMatrix) {
  const int numFeatures = 5;
  const int numRows = 5000;  // more than the first matrix of 1024 rows
  std::vector<double> features(static_cast<std::size_t>(numRows) * numFeatures);
  for (std::size_t i = 0; i < features.size(); ++i) {
    features[i] = static_cast<double>(i);
  }
  AlgIn data(numFeatures + 1);
  for (int i = 0; i < numRows; ++i) {
    data.add(&features[static_cast<std::size_t>(i) * numFeatures], i % 2 ? 1 : -1);
  }
  AlgIn copy(data);
  AlgIn assigned(numFeatures + 1);
  assigned = data;
  AlgIn moved(std::move(copy));
  for (const AlgIn* set : {&data, &assigned, &moved}) {
    ASSERT_EQ(numRows, set->m);
    for (int i = 0; i < numRows; ++i) {
      for (int k = 0; k < numFeatures; ++k) {
        ASSERT_EQ(features[static_cast<std::size_t>(i) * numFeatures + k],
                  set->row(i)[k]);
      }
      ASSERT_EQ(1.0, set->row(i)[numFeatures]);
    }
  }
  EXPECT_NE(data.row(0), assigned.row(0));
}

// The solver gives the same weights and outputs on copied rows and on rows
// in place.
TEST(TrainingMemoryTest, inPlaceRowsGiveTheSameSolution) {
  const int numFeatures = 6;
  const int n = numFeatures + 1;
  const int m = 900;
  std::vector<double> storage;
  AlgIn inPlace(n, true);
  makeTrainingSet(m, numFeatures, 99u, storage, inPlace);
  AlgIn copied(n);
  for (int i = 0; i < m; ++i) {
    copied.add(inPlace.row(i), inPlace.Y[i]);
  }
  copied.negatives = inPlace.negatives;
  copied.positives = inPlace.positives;
  for (bool directSolve : {true, false}) {
    options opts;
    opts.lambda = 1.0;
    opts.lambda_u = 1.0;
    opts.epsilon = EPSILON;
    opts.cgitermax = CGITERMAX;
    opts.mfnitermax = MFNITERMAX;
    opts.directSolve = directSolve;
    std::vector<double> w1(n, 0.0), o1(m, 0.0), w2(n, 0.0), o2(m, 0.0);
    L2_SVM_MFN(copied, opts, w1, o1, 1.0, 2.0);
    L2_SVM_MFN(inPlace, opts, w2, o2, 1.0, 2.0);
    EXPECT_NE(0.0, w1[0]);
    for (int k = 0; k < n; ++k) {
      EXPECT_EQ(w1[k], w2[k]);
    }
    for (int i = 0; i < m; ++i) {
      EXPECT_EQ(o1[i], o2[i]);
    }
  }
}
