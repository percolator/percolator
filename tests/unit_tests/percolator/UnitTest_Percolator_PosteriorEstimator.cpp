/*******************************************************************************
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
 * Unit tests for reading the input files of qvality and for its q values
 * from p values.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "Globals.h"
#include "LogisticRegression.h"
#include "MyException.h"
#include "PosteriorEstimator.h"
#include "TmpDir.h"

namespace {

class ScoreReader : public PosteriorEstimator {
 public:
  using PosteriorEstimator::readScores;
};

class Estimator : public PosteriorEstimator {
 public:
  using PosteriorEstimator::estimate;
};

// the knots of a P-IRLS fit after LogisticRegression::setData
class Knots : public LogisticRegression {
 public:
  std::size_t size() const { return x.size(); }
};

class KnotMerger : public LogisticRegression {
 public:
  using LogisticRegression::countNearTies;
  using LogisticRegression::mergeEqualKnots;
};

// the knots that are left after mergeEqualKnots
std::vector<double> mergedKnots(std::vector<double> x) {
  std::vector<double> y(x.size(), 1.0), m(x.size(), 2.0);
  KnotMerger::mergeEqualKnots(x, y, m);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_EQ(x.size(), m.size());
  return x;
}

// uniform in [0, 1), from 53 bits of the generator
double uniform53(std::mt19937_64& rng) {
  return static_cast<double>(rng() >> 11) * 0x1p-53;
}

std::string writeFile(const TmpDir& dir, const std::string& text) {
  std::filesystem::path path =
      std::filesystem::path(dir.getPath()) / "scores.txt";
  std::ofstream(path.string()) << text;
  return path.string();
}

// Sets the verbosity and the no-terminate flag, and restores them at the end
// of the test.
class GlobalFlags {
 public:
  GlobalFlags(int verbose, bool noTerminate)
      : verbose_(Globals::getInstance()->getVerbose()),
        noTerminate_(Globals::getInstance()->getNoTerminate()) {
    Globals::getInstance()->setVerbose(verbose);
    Globals::getInstance()->setNoTerminate(noTerminate);
  }
  ~GlobalFlags() {
    Globals::getInstance()->setVerbose(verbose_);
    Globals::getInstance()->setNoTerminate(noTerminate_);
  }

 private:
  int verbose_;
  bool noTerminate_;
};

}  // namespace

TEST(PosteriorEstimatorTest, ReadScores) {
  std::vector<double> scores;
  TmpDir tmpDir;
  std::string fn = writeFile(tmpDir, "1.5 -2\n3e2\t0.25\n\n");
  ScoreReader::readScores(fn, "target score", scores);
  EXPECT_EQ((std::vector<double>{1.5, -2.0, 300.0, 0.25}), scores);

  // no newline at the end
  scores.clear();
  fn = writeFile(tmpDir, "1\n2");
  ScoreReader::readScores(fn, "target score", scores);
  EXPECT_EQ((std::vector<double>{1.0, 2.0}), scores);

  // reading used to stop silently at the first token that is not a number
  for (const char* text : {"1\n2\nnan\n3\n", "1\nabc\n", "1\n2.5x\n",
                           "inf\n", "1e999\n", "", " \n\n"}) {
    scores.clear();
    fn = writeFile(tmpDir, text);
    EXPECT_THROW(ScoreReader::readScores(fn, "target score", scores),
                 MyException) << text;
  }
  std::filesystem::remove(fn);
  scores.clear();
  EXPECT_THROW(ScoreReader::readScores(fn, "null score", scores), MyException);
}

// The help, also the one printed for a wrong number of arguments, used to
// end the process inside the option parser. parseOptions returns true and
// run does nothing, so qvality ends with status 0 as before.
// qvality -d on a p value file prints the p values and as many generated
// decoys; the q values used to be those of the p values only, and the
// output read past their end.
TEST(PosteriorEstimatorTest, QValuesFromPWithDecoys) {
  const std::vector<double> p = {0.01, 0.02, 0.5, 0.9};
  std::vector<double> targetQ;
  PosteriorEstimator::getQValuesFromP(1.0, p, targetQ);
  ASSERT_EQ(p.size(), targetQ.size());
  // the decoys of PosteriorEstimator::run, (2i + 1) / (2m), merged in
  const std::vector<std::pair<double, bool> > combined = {
      {0.01, true},   {0.02, true},  {0.125, false}, {0.375, false},
      {0.5, true},    {0.625, false}, {0.875, false}, {0.9, true}};
  const std::vector<double> q =
      PosteriorEstimator::getAllQValuesFromP(1.0, p, targetQ, combined);
  // m * d / #{p <= d}, limited by the next target: 4 * 0.125 / 2,
  // min(4 * 0.375 / 2, 2/3), 4 * 0.625 / 3, min(4 * 0.875 / 3, 0.9)
  const std::vector<double> expected = {0.04, 0.04, 0.25,      2.0 / 3.0,
                                        2.0 / 3.0, 2.5 / 3.0, 0.9, 0.9};
  ASSERT_EQ(expected.size(), q.size());
  for (std::size_t i = 0; i < q.size(); ++i) {
    EXPECT_NEAR(expected[i], q[i], 1e-12) << "entry " << i;
    if (i > 0) {
      EXPECT_LE(q[i - 1], q[i]) << "entry " << i;
    }
  }

  // a decoy before the first p value, and pi0
  const std::vector<double> late = {0.9, 0.95};
  targetQ.clear();
  PosteriorEstimator::getQValuesFromP(0.5, late, targetQ);
  const std::vector<double> q2 = PosteriorEstimator::getAllQValuesFromP(
      0.5, late, targetQ,
      {{0.25, false}, {0.75, false}, {0.9, true}, {0.95, true}});
  ASSERT_EQ(4u, q2.size());
  EXPECT_NEAR(0.25, q2[0], 1e-12);
  EXPECT_NEAR(0.475, q2[1], 1e-12);  // limited by the targets, 0.5 * 0.95
  EXPECT_NEAR(0.475, q2[2], 1e-12);
  EXPECT_NEAR(0.475, q2[3], 1e-12);
}

// The P-IRLS PEP fit (protein PEPs, --irls-pep) on fewer than 500 entries,
// one per bin, with two scores 3e-6 apart. In double precision the error of
// the linear solve in each Newton step was as large as the tolerance on the
// step size (1e-8), so the number of steps and the PEPs depended on
// rounding: moving one score by 1 ulp changed the PEPs by 2.3e-7 relative.
// With the step in double-double they do not change.
TEST(PosteriorEstimatorTest, IrlsPepStableUnderOneUlpChange) {
  PosteriorEstimator::resetToDefaults();
  std::mt19937_64 rng(7);
  auto uniform = [&rng]() { return static_cast<double>(rng() >> 11) * 0x1p-53; };
  const int n = 240;
  std::vector<std::pair<double, bool> > combined;
  for (int i = 0; i < n; ++i) {
    combined.emplace_back(-2.0 + 4.0 * uniform(), true);
  }
  combined.emplace_back(combined[7].first + 3e-6, true);
  std::sort(combined.begin(), combined.end());
  // lower scores are better; decoys become more frequent towards the top
  for (auto& entry : combined) {
    double pDecoy = 0.45 * std::min(1.0, std::max(0.0, (entry.first + 1.0) / 3.0));
    entry.second = !(uniform() < pDecoy);
  }
  std::vector<std::pair<double, bool> > moved = combined;
  const std::size_t ix = combined.size() / 2;
  moved[ix].first = std::nextafter(moved[ix].first, 10.0);
  ASSERT_LT(moved[ix].first, moved[ix + 1].first);

  std::vector<double> peps, pepsMoved;
  PosteriorEstimator::estimatePEP(combined, false, 1.0, peps, true, 1.0, true);
  PosteriorEstimator::estimatePEP(moved, false, 1.0, pepsMoved, true, 1.0,
                                  true);
  ASSERT_EQ(combined.size(), peps.size());
  ASSERT_EQ(peps.size(), pepsMoved.size());
  double maxRel = 0.0;
  for (std::size_t i = 0; i < peps.size(); ++i) {
    ASSERT_GT(peps[i], 0.0) << i;
    ASSERT_LE(peps[i], 1.0) << i;
    // also fails for NaN, which std::max below would pass over
    ASSERT_GT(pepsMoved[i], 0.0) << i;
    ASSERT_LE(pepsMoved[i], 1.0) << i;
    if (i > 0) {
      EXPECT_LE(peps[i - 1], peps[i]) << i;
    }
    maxRel = std::max(maxRel, std::fabs(peps[i] - pepsMoved[i]) / peps[i]);
  }
  EXPECT_LT(maxRel, 1e-10);
}

// Two equal knots give an interval of length 0, and the step of the P-IRLS
// fit was not a number. The bound on g turned it into -35, the next step was
// 0 and the fit was taken as converged: every PEP was exp(-35) = 6.3e-16,
// without a warning. Now the step fails, g keeps its start values, and the
// warning is printed.
TEST(PosteriorEstimatorTest, IrlsFitFailsAtEqualKnots) {
  GlobalFlags flags(1, false);
  const std::vector<double> x = {-1.0, 0.0, 0.0, 1.0};
  const std::vector<double> y = {2.0, 1.0, 1.0, 0.0};
  const std::vector<double> m = {2.0, 2.0, 2.0, 2.0};
  LogisticRegression lr;
  lr.setData(x, y, m);
  testing::internal::CaptureStderr();
  lr.roughnessPenaltyIRLS();
  const std::string warning = testing::internal::GetCapturedStderr();
  EXPECT_EQ("Warning: the P-IRLS fit of the PEPs did not converge (4 bins); "
            "the PEPs may be inaccurate.\n",
            warning);
  EXPECT_FALSE(lr.iterativeReweightedLeastSquares(1.0));
  std::vector<double> g;
  lr.predict(x, g);
  ASSERT_EQ(x.size(), g.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    // the start value of LogisticRegression::initg
    EXPECT_NEAR(std::log((y[i] + 0.05) / (m[i] - y[i] + 0.05)), g[i], 1e-12)
        << i;
  }
}

// With --no-terminate, fewer than 4 bins are padded to 4. The padded bins
// had the score of the last bin (the offset of 1e-20 was lost unless the
// score was close to 0), the fit failed as above, and every PEP was 6.3e-16,
// also for a list of decoys only. Tied scores, a single bin, were given
// noise of the same size, with the same result. Now the padded bins have
// distinct scores; a single bin gives every score the same PEP, and the
// scores are not changed.
TEST(PosteriorEstimatorTest, NoTerminatePadsWithDistinctScores) {
  GlobalFlags flags(0, true);
  // sorted by increasing score, the best first, as protein scores
  std::vector<std::pair<double, bool> > decoys(12, std::make_pair(1.0, false));
  std::vector<double> peps;
  PosteriorEstimator::estimatePEP(decoys, false, 1.0, peps, true, 1.0, true);
  ASSERT_EQ(12u, peps.size());
  for (std::size_t i = 0; i < peps.size(); ++i) {
    EXPECT_EQ(1.0, peps[i]) << i;
    EXPECT_EQ(1.0, decoys[i].first) << i;
  }

  // 3 decoys and 9 targets with the same score: the odds 3 / 9
  std::vector<std::pair<double, bool> > tied(3, std::make_pair(2.5, false));
  tied.resize(12, std::make_pair(2.5, true));
  PosteriorEstimator::estimatePEP(tied, false, 1.0, peps, true, 1.0, true);
  ASSERT_EQ(12u, peps.size());
  for (std::size_t i = 0; i < peps.size(); ++i) {
    EXPECT_NEAR(1.0 / 3.0, peps[i], 1e-6) << i;
  }

  // three bins, two targets and a decoy with a worse score
  std::vector<std::pair<double, bool> > three = {
      {1.0, true}, {2.0, true}, {3.0, false}};
  PosteriorEstimator::estimatePEP(three, false, 1.0, peps, true, 1.0, true);
  ASSERT_EQ(3u, peps.size());
  EXPECT_GT(peps[1], peps[0]);
  EXPECT_LT(peps[1], 1e-6);
  EXPECT_EQ(1.0, peps[2]);
}

// Equal knots are merged, with their y and m summed. Knots that are only
// close together are not, and nothing is merged if fewer than 4 knots would
// remain.
TEST(PosteriorEstimatorTest, MergeEqualKnots) {
  std::vector<double> x = {-1.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  std::vector<double> m = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
  KnotMerger::mergeEqualKnots(x, y, m);
  EXPECT_EQ((std::vector<double>{-1.0, 0.0, 1.0, 2.0}), x);
  EXPECT_EQ((std::vector<double>{1.0, 5.0, 4.0, 18.0}), y);
  EXPECT_EQ((std::vector<double>{10.0, 50.0, 40.0, 180.0}), m);

  const std::vector<double> oneUlp = {-1.0, 0.0, 1.0,
                                      std::nextafter(1.0, 2.0), 2.0};
  EXPECT_EQ(oneUlp, mergedKnots(oneUlp));
  const std::vector<double> distinct = {-1.0, 0.0, 1.0, 2.0};
  EXPECT_EQ(distinct, mergedKnots(distinct));
  const std::vector<double> tiedOfFour = {0.0, 1.0, 1.0, 2.0};
  EXPECT_EQ(tiedOfFour, mergedKnots(tiedOfFour));
  const std::vector<double> tiedOfFive = {0.0, 1.0, 1.0, 1.0, 2.0};
  EXPECT_EQ(tiedOfFive, mergedKnots(tiedOfFive));
}

// Two scores near 1e9 one ulp apart are two bins, and the log transform of
// BaseSpline::setData gives them the same knot. The fit could take no step
// and ended with the warning and the start values of g. The knots are now
// merged after the transform: the fit converges, and the PEPs are those of
// the same entries with the two scores equal, which share one bin.
TEST(PosteriorEstimatorTest, IrlsPepWithEqualKnotsAfterTheLog) {
  PosteriorEstimator::resetToDefaults();
  GlobalFlags flags(1, false);
  std::mt19937_64 rng(3);
  // lower scores are better; decoys become more frequent towards the top
  std::vector<std::pair<double, bool> > combined;
  for (int i = 0; i < 400; ++i) {
    combined.emplace_back(1e9 + 1000.0 * (i + 0.4 * uniform53(rng)),
                          !(uniform53(rng) < 0.45 * i / 400.0));
  }
  // a score whose next double has the same log, on this platform
  double s = combined[200].first;
  while (std::log(s) != std::log(std::nextafter(s, 2e9))) {
    s = std::nextafter(s, 2e9);
  }
  combined[200].first = s;
  std::vector<std::pair<double, bool> > tied = combined;
  combined[201].first = std::nextafter(s, 2e9);
  tied[201].first = s;
  ASSERT_LT(combined[201].first, combined[202].first);

  Knots knots, tiedKnots;
  Estimator::estimate(combined, knots, false, 1.0, true);
  Estimator::estimate(tied, tiedKnots, false, 1.0, true);
  EXPECT_EQ(399u, knots.size());
  EXPECT_EQ(399u, tiedKnots.size());

  std::vector<double> peps, pepsTied;
  testing::internal::CaptureStderr();
  PosteriorEstimator::estimatePEP(combined, false, 1.0, peps, true, 1.0, true);
  PosteriorEstimator::estimatePEP(tied, false, 1.0, pepsTied, true, 1.0,
                                  true);
  const std::string warning = testing::internal::GetCapturedStderr();
  EXPECT_EQ(std::string::npos, warning.find("did not converge")) << warning;
  ASSERT_EQ(combined.size(), peps.size());
  EXPECT_EQ(pepsTied, peps);
}

// Pairs of knots less than 1e-9 times the median of the gaps around them
// apart are counted as nearly equal.
TEST(PosteriorEstimatorTest, CountNearTies) {
  // -2 to 2, 0.01 apart
  std::vector<double> x;
  for (int i = 0; i <= 400; ++i) {
    x.push_back((i - 200) / 100.0);
  }
  EXPECT_EQ(0u, KnotMerger::countNearTies(x));
  for (double offset : {1e-12, 1e-10}) {
    std::vector<double> pair = x;
    pair.insert(pair.begin() + 301, pair[300] + offset);
    EXPECT_EQ(offset < 1e-11 ? 1u : 0u, KnotMerger::countNearTies(pair))
        << offset;
  }
  // 0 and 1e-20
  std::vector<double> nearZero = x;
  nearZero.insert(nearZero.begin() + 201, 1e-20);
  EXPECT_EQ(1u, KnotMerger::countNearTies(nearZero));
  // 24 knots 1 ulp apart from 1.5 count at both ends of the run, also
  // where they are more than half of the knots
  std::vector<double> run(x.begin(), x.begin() + 351);
  for (int j = 1; j < 24; ++j) {
    run.push_back(std::nextafter(run.back(), 2.0));
  }
  run.insert(run.end(), x.begin() + 351, x.end());
  EXPECT_EQ(2u, KnotMerger::countNearTies(run));
  std::vector<double> mostly(run.begin() + 340, run.begin() + 384);
  ASSERT_EQ(1.4, mostly.front());
  EXPECT_EQ(2u, KnotMerger::countNearTies(mostly));
  // all knots in runs 1 ulp apart from 1, 1.04, ..., 4.96: in runs of 2
  // every pair counts; in runs of 3 or 4 the gaps around each pair are small
  // too, and none does
  for (int length : {2, 3, 4}) {
    std::vector<double> runs;
    for (int i = 0; i < 100; ++i) {
      runs.push_back(1.0 + i / 25.0);
      for (int j = 1; j < length; ++j) {
        runs.push_back(std::nextafter(runs.back(), 10.0));
      }
    }
    EXPECT_EQ(length == 2 ? 100u : 0u, KnotMerger::countNearTies(runs))
        << length;
  }
}

// Two of 240 scores 1e-12 apart (5.6e-11 times the median gap around them)
// are nearly equal, and the warning is printed once, before the fit; 1e-8
// apart they are not. Only this warning is checked: with the two scores
// 1e-13 apart, the fit converged or reached its limit depending on the
// rounding of exp and log.
TEST(PosteriorEstimatorTest, IrlsPepWarnsOfNearlyEqualScores) {
  PosteriorEstimator::resetToDefaults();
  GlobalFlags flags(1, false);
  std::mt19937_64 rng(7);
  // lower scores are better; decoys become more frequent towards the top
  std::vector<std::pair<double, bool> > combined;
  for (int i = 0; i < 240; ++i) {
    combined.emplace_back(-2.0 + 4.0 * uniform53(rng), true);
  }
  std::sort(combined.begin(), combined.end());
  for (auto& entry : combined) {
    const double pDecoy =
        0.45 * std::min(1.0, std::max(0.0, (entry.first + 1.0) / 3.0));
    entry.second = !(uniform53(rng) < pDecoy);
  }
  const std::string warning =
      "Warning: 1 pairs of scores in the PEP fit are nearly equal; the "
      "PEPs may be inaccurate.\n";
  for (double offset : {1e-12, 1e-8}) {
    std::vector<std::pair<double, bool> > entries = combined;
    entries[121].first = entries[120].first + offset;
    ASSERT_LT(entries[121].first, entries[122].first);
    std::vector<double> peps;
    testing::internal::CaptureStderr();
    PosteriorEstimator::estimatePEP(entries, false, 1.0, peps, true, 1.0,
                                    true);
    const std::string err = testing::internal::GetCapturedStderr();
    if (offset < 1e-10) {
      EXPECT_EQ(0u, err.find(warning)) << err;
      EXPECT_EQ(std::string::npos, err.find("nearly equal", warning.size()))
          << err;
    } else {
      EXPECT_EQ(std::string::npos, err.find("nearly equal")) << err;
    }
  }
}

// Probabilities with ties at 1.0 and one score 1e-12 below. The logit
// transform is chosen from all scores, before equal knots are merged, and
// has an offset at 1, so that 1.0 has a finite knot: no PEP is 0 or above
// 1, and fewer than half of them are 1.
TEST(PosteriorEstimatorTest, IrlsPepWithTiesAtOne) {
  PosteriorEstimator::resetToDefaults();
  std::mt19937_64 rng(5);
  std::vector<std::pair<double, bool> > combined;
  for (int i = 0; i < 190; ++i) {
    const double u = uniform53(rng);
    combined.emplace_back(1.0 - 0.9 * u * u * u, true);
  }
  combined[0].first = 1.0 - 1e-12;
  combined.resize(200, std::make_pair(1.0, true));
  for (int i = 0; i < 200; ++i) {
    combined.emplace_back(0.95 * uniform53(rng), false);
  }
  for (int decoysAtOne : {0, 6}) {
    std::vector<std::pair<double, bool> > entries = combined;
    entries.resize(400 + decoysAtOne, std::make_pair(1.0, false));
    // higher scores are better
    std::sort(entries.begin(), entries.end(),
              std::greater<std::pair<double, bool> >());
    std::vector<double> peps;
    PosteriorEstimator::estimatePEP(entries, false, 1.0, peps, false, 1.0,
                                    false);
    ASSERT_EQ(200u, peps.size());
    int outside = 0, ones = 0;
    for (double pep : peps) {
      outside += !(pep > 0.0 && pep <= 1.0);
      ones += pep == 1.0;
    }
    EXPECT_EQ(0, outside) << decoysAtOne;
    EXPECT_LT(ones, 100) << decoysAtOne;
    EXPECT_LT(peps[0], decoysAtOne == 0 ? 1e-6 : 0.01) << decoysAtOne;
  }
}

TEST(PosteriorEstimatorTest, HelpReturnsWithoutEndingTheProcess) {
  const std::vector<std::vector<std::string> > commands = {
      {"qvality", "--help", "--no-such-option"},
      {"qvality"},
      {"qvality", "a.txt", "b.txt", "c.txt"}};
  for (std::vector<std::string> args : commands) {
    std::vector<char*> argv;
    for (std::string& arg : args) argv.push_back(&arg[0]);
    PosteriorEstimator estimator;
    testing::internal::CaptureStderr();
    bool parsed = false;
    EXPECT_NO_THROW(parsed = estimator.parseOptions(
                        static_cast<int>(args.size()), argv.data()));
    int ran = 0;
    EXPECT_NO_THROW(ran = estimator.run());
    const std::string help = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(parsed) << args.size();
    EXPECT_TRUE(ran) << args.size();
    EXPECT_NE(std::string::npos, help.find("Options:")) << args.size();
  }
}
