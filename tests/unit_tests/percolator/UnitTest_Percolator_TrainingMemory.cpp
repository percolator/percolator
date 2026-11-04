/*
 * Unit tests for the memory reporting of the stage lines and for the
 * memory handling of the SVM training.
 */

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

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
