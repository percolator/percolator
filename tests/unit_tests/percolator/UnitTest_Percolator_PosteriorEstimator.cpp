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

#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "MyException.h"
#include "PosteriorEstimator.h"
#include "TmpDir.h"

namespace {

class ScoreReader : public PosteriorEstimator {
 public:
  using PosteriorEstimator::readScores;
};

std::string writeFile(const TmpDir& dir, const std::string& text) {
  std::filesystem::path path =
      std::filesystem::path(dir.getPath()) / "scores.txt";
  std::ofstream(path.string()) << text;
  return path.string();
}

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
    if (i > 0) EXPECT_LE(q[i - 1], q[i]) << "entry " << i;
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
