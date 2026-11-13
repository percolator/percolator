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

// A training set of more PSMs than the 100,000 from which the rows are
// counted and copied in one chunk per thread equals the one made by
// appending the selected PSMs one by one, with and without the best
// positive per spectrum, also when it is built again in the same AlgIn.
TEST(TrainingMemoryTest, largeTrainingSetEqualsRowByRowAppend) {
  const unsigned int numFeatures = 5;
  FeatureNames::resetNumFeatures();
  FeatureNames::setNumFeatures(numFeatures);
  const std::size_t numPsms = 250007;
  std::vector<double> features(numPsms * numFeatures);
  std::vector<PSMDescription> psms(numPsms);
  Scores scores(false);
  unsigned int seed = 20261u;
  for (std::size_t i = 0; i < numPsms; ++i) {
    for (unsigned int k = 0; k < numFeatures; ++k) {
      seed = seed * 1103515245u + 12345u;
      features[i * numFeatures + k] = static_cast<double>(seed >> 8) / 1024.0;
    }
    psms[i].features = &features[i * numFeatures];
    psms[i].scan = static_cast<unsigned int>(i / 3);  // three PSMs per scan
    psms[i].setLineNr(static_cast<unsigned int>(i));
    seed = seed * 1103515245u + 12345u;
    const bool isDecoy = (seed >> 16) % 3 == 0;
    ScoreHolder sh(static_cast<double>(numPsms - i),
                   isDecoy ? LabelType::DECOY : LabelType::TARGET, &psms[i]);
    seed = seed * 1103515245u + 12345u;
    sh.q = (seed >> 16) % 4 == 0 ? 0.5 : 0.001;
    scores.addScoreHolder(sh);
  }
  for (bool bestPositive : {false, true}) {
    // the reference: the selected PSMs, one by one
    std::vector<const PSMDescription*> expected;
    for (const ScoreHolder& sh : scores) {
      if (sh.isDecoy()) expected.push_back(sh.pPSM);
    }
    const std::size_t numNegatives = expected.size();
    unsigned int lastScan = static_cast<unsigned int>(-1);
    for (const ScoreHolder& sh : scores) {
      if (!sh.isTarget() || sh.q > 0.01) continue;
      // the scores decrease, so the first target of a scan is its best
      if (bestPositive && sh.pPSM->scan == lastScan) continue;
      lastScan = sh.pPSM->scan;
      expected.push_back(sh.pPSM);
    }
    ASSERT_GT(numNegatives, 50000u);
    ASSERT_GT(expected.size() - numNegatives, 50000u);

    AlgIn data(static_cast<int>(numFeatures) + 1);
    const double* firstRow = nullptr;
    for (int round = 0; round < 2; ++round) {
      scores.generateNegativeTrainingSet(data, 1.0);
      ASSERT_EQ(static_cast<int>(numNegatives), data.m);
      scores.generatePositiveTrainingSet(data, 0.01, 1.0, bestPositive);
      ASSERT_EQ(static_cast<int>(expected.size()), data.m);
      ASSERT_EQ(static_cast<int>(numNegatives), data.negatives);
      ASSERT_EQ(static_cast<int>(expected.size() - numNegatives),
                data.positives);
      ASSERT_EQ(expected.size(), data.Y.size());
      // the matrix is allocated once
      if (round == 0) firstRow = data.row(0);
      ASSERT_EQ(firstRow, data.row(0));
      std::size_t numWrong = 0;
      for (std::size_t i = 0; i < expected.size(); ++i) {
        const double* row = data.row(static_cast<int>(i));
        for (unsigned int k = 0; k < numFeatures; ++k) {
          if (row[k] != expected[i]->features[k]) ++numWrong;
        }
        if (row[numFeatures] != 1.0) ++numWrong;
        if (data.Y[i] != (i < numNegatives ? -1.0 : 1.0)) ++numWrong;
      }
      ASSERT_EQ(0u, numWrong);
    }
  }
  FeatureNames::resetNumFeatures();
}

// reserveTrainingSet makes room for the decoys and the accepted targets
// and an eighth more, not for all PSMs, and keeps a buffer that is large
// enough.
TEST(TrainingMemoryTest, reservesRoomForTheTrainingSetOnly) {
  const unsigned int numFeatures = 3;
  FeatureNames::resetNumFeatures();
  FeatureNames::setNumFeatures(numFeatures);
  for (bool inPlace : {false, true}) {
    const unsigned int numPsms = 4000;
    std::vector<double> features(numPsms * (numFeatures + 1), 1.0);
    std::vector<PSMDescription> psms(numPsms);
    Scores scores(false);
    for (unsigned int i = 0; i < numPsms; ++i) {
      psms[i].features = &features[i * (numFeatures + 1)];
      psms[i].scan = i;
      // 1000 decoys; 3000 targets, of which 600 with q <= 0.01
      ScoreHolder sh(1.0 * (numPsms - i),
                     i % 4 == 3 ? LabelType::DECOY : LabelType::TARGET,
                     &psms[i]);
      sh.q = i < 800 ? 0.01 : 0.5;
      scores.addScoreHolder(sh);
    }
    AlgIn data(static_cast<int>(numFeatures) + 1, inPlace);
    scores.reserveTrainingSet(data, 0.01);
    EXPECT_EQ(0, data.m);
    EXPECT_GE(data.capacity(), 1600u + 200u);
    EXPECT_LT(data.capacity(), 2000u);
    const std::size_t capacity = data.capacity();
    scores.generateNegativeTrainingSet(data, 1.0);
    scores.generatePositiveTrainingSet(data, 0.01, 1.0, false);
    EXPECT_EQ(1000, data.negatives);
    EXPECT_EQ(600, data.positives);
    const double* firstRow = data.row(0);
    // the next, smaller, training set uses the same buffer
    scores.reserveTrainingSet(data, 0.001);
    EXPECT_EQ(0, data.m);
    EXPECT_EQ(capacity, data.capacity());
    scores.generateNegativeTrainingSet(data, 1.0);
    scores.generatePositiveTrainingSet(data, 0.001, 1.0, false);
    EXPECT_EQ(1000, data.m);
    EXPECT_EQ(firstRow, data.row(0));
    // a larger one gets a new buffer
    scores.reserveTrainingSet(data, 0.5);
    EXPECT_GE(data.capacity(), 4000u + 500u);
    scores.generateNegativeTrainingSet(data, 1.0);
    scores.generatePositiveTrainingSet(data, 0.5, 1.0, false);
    EXPECT_EQ(4000, data.m);
  }
  FeatureNames::resetNumFeatures();
}

// resizeExamples and set give the examples of add, clear keeps the matrix,
// and examples written before a resize survive it.
TEST(TrainingMemoryTest, examplesSetInPlaceEqualAppendedOnes) {
  const int numFeatures = 4;
  const int numRows = 3000;
  std::vector<double> features(static_cast<std::size_t>(numRows) * numFeatures);
  for (std::size_t i = 0; i < features.size(); ++i) {
    features[i] = 0.25 * static_cast<double>(i);
  }
  AlgIn appended(numFeatures + 1);
  for (int i = 0; i < numRows; ++i) {
    appended.add(&features[static_cast<std::size_t>(i) * numFeatures],
                 i % 3 ? 1 : -1);
  }
  AlgIn data(numFeatures + 1);
  for (int round = 0; round < 2; ++round) {
    // in two parts, the second one written backwards
    data.resizeExamples(1000);
    ASSERT_EQ(1000, data.m);
    for (std::size_t i = 0; i < 1000; ++i) {
      data.set(i, &features[i * numFeatures], i % 3 ? 1 : -1);
    }
    const double* firstRow = data.row(0);
    data.resizeExamples(static_cast<std::size_t>(numRows));
    if (round == 1) {
      EXPECT_EQ(firstRow, data.row(0));
    }
    for (std::size_t i = static_cast<std::size_t>(numRows); i-- > 1000;) {
      data.set(i, &features[i * numFeatures], i % 3 ? 1 : -1);
    }
    ASSERT_EQ(appended.m, data.m);
    ASSERT_EQ(appended.Y, data.Y);
    for (int i = 0; i < numRows; ++i) {
      for (int k = 0; k <= numFeatures; ++k) {
        ASSERT_EQ(appended.row(i)[k], data.row(i)[k]);
      }
    }
    data.positives = 1;
    data.negatives = 2;
    data.clear();
    EXPECT_EQ(0, data.m);
    EXPECT_EQ(0, data.positives);
    EXPECT_EQ(0, data.negatives);
    EXPECT_TRUE(data.Y.empty());
  }
  // in place, the examples are the given rows
  std::vector<double> rows(3 * (numFeatures + 1), 1.0);
  AlgIn inPlace(numFeatures + 1, true);
  inPlace.resizeExamples(3);
  for (std::size_t i = 3; i-- > 0;) {
    inPlace.set(i, &rows[i * (numFeatures + 1)], -1);
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(&rows[static_cast<std::size_t>(i) * (numFeatures + 1)],
              inPlace.row(i));
    EXPECT_EQ(-1.0, inPlace.Y[static_cast<std::size_t>(i)]);
  }
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
  // a set that was emptied and one that was moved from can be copied,
  // assigned to and filled again
  data.clear();
  AlgIn emptyCopy(data);
  EXPECT_EQ(0, emptyCopy.m);
  emptyCopy.add(&features[0], 1);
  EXPECT_EQ(features[1], emptyCopy.row(0)[1]);
  EXPECT_EQ(0, copy.m);
  copy.add(&features[numFeatures], -1);
  EXPECT_EQ(features[numFeatures], copy.row(0)[0]);
  copy = assigned;
  ASSERT_EQ(numRows, copy.m);
  EXPECT_EQ(features[features.size() - 1],
            copy.row(numRows - 1)[numFeatures - 1]);
  EXPECT_EQ(numRows, assigned.m);
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

// The normal equations of the direct solve, which are built from four rows
// per pass, equal bit for bit those built from one row per pass, in which
// every element receives the terms of the rows in the same order: for
// numbers of active rows with every remainder of four, on copied rows and
// on rows in place.
TEST(TrainingMemoryTest, normalEquationsOfFourRowsPerPassAreBitwiseEqual) {
  const int numFeatures = 9;
  const int n = numFeatures + 1;
  const int m = 61;
  std::vector<double> storage;
  AlgIn inPlace(n, true);
  makeTrainingSet(m, numFeatures, 4242u, storage, inPlace);
  // values with full mantissas, so that the order of the additions shows
  for (int i = 0; i < m; ++i) {
    for (int k = 0; k < numFeatures; ++k) {
      double& value = storage[static_cast<std::size_t>(i) * n + k];
      value = value / 3.0 + 1e-9 * (i * 31 + k);
    }
  }
  AlgIn copied(n);
  for (int i = 0; i < m; ++i) {
    copied.add(inPlace.row(i), inPlace.Y[i]);
  }
  const double cpos = 0.7, cneg = 2.1;
  for (const AlgIn* data : {&copied, &inPlace}) {
    for (int active = 0; active <= 23; ++active) {
      // active rows that skip through the set, in an order of their own
      std::vector<int> J(static_cast<std::size_t>(active));
      for (int i = 0; i < active; ++i) J[i] = (i * 37 + 5) % m;
      std::vector<double> expectedG(static_cast<std::size_t>(n) * n, 0.0);
      std::vector<double> expectedB(static_cast<std::size_t>(n), 0.0);
      for (int i = 0; i < active; ++i) {
        const double* x = data->row(J[i]);
        const double y = data->Y[J[i]];
        const double c = y == 1 ? cpos : cneg;
        for (int j = 0; j < n; ++j) {
          const double cx = c * x[j];
          expectedB[j] += cx * y;
          for (int k = j; k < n; ++k) {
            expectedG[static_cast<std::size_t>(j) * n + k] += cx * x[k];
          }
        }
      }
      // G and b are overwritten
      std::vector<double> G(static_cast<std::size_t>(n) * n, -7.0);
      std::vector<double> b(static_cast<std::size_t>(n), -7.0);
      accumulateNormalEquations(*data, J.data(), active, cpos, cneg, G.data(),
                                b.data());
      for (int j = 0; j < n; ++j) {
        ASSERT_EQ(expectedB[j], b[j]) << active << " rows, b " << j;
        for (int k = j; k < n; ++k) {
          ASSERT_EQ(expectedG[static_cast<std::size_t>(j) * n + k],
                    G[static_cast<std::size_t>(j) * n + k])
              << active << " rows, G " << j << " " << k;
        }
      }
    }
  }
}
