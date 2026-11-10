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
 * Unit tests for reading the input files of qvality.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>
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
