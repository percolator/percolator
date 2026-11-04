/*
 * Unit tests for the memory reporting of the stage lines and for the
 * memory handling of the SVM training.
 */

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

#include "FeatureNames.h"
#include "PSMDescription.h"
#include "Scores.h"
#include "Timer.h"

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
