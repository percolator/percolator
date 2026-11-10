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
 * Runs Percolator several times in one process, as Crux does, and checks
 * that a run gives the same result after another run as on its own.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "Caller.h"
#include "TmpDir.h"
#include "Globals.h"

namespace {

namespace fs = std::filesystem;

// a concatenated search of 2,000 spectra with one PSM each; two thirds of
// the targets score higher in two of the four features, so that the SVM
// improves on the best single feature
void writeInput(const std::string& fileName) {
  std::mt19937 rng(17);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::ofstream out(fileName);
  out << "SpecId\tLabel\tScanNr\tExpMass\tf1\tf2\tf3\tf4\tPeptide\tProteins\n";
  const char* residues = "ACDEFGHIKLMNPQRSTVWY";
  for (int i = 0; i < 2000; ++i) {
    const bool isTarget = rng() % 2 == 0;
    const double shift = isTarget && rng() % 3 != 0 ? 2.0 : 0.0;
    std::string peptide = "K.";
    for (int k = i + 10000; k > 0; k /= 20) peptide += residues[k % 20];
    peptide += "K.A";  // tryptic, as the protein inference expects
    out << (isTarget ? "target_" : "decoy_") << i << '\t'
        << (isTarget ? 1 : -1) << '\t' << i + 1 << '\t'
        << 1000.0 + i * 0.5 << '\t' << noise(rng) + shift << '\t'
        << noise(rng) + shift << '\t' << noise(rng) << '\t'
        << noise(rng) << '\t' << peptide << '\t'
        << (isTarget ? "prot_" : "decoy_prot_") << i % 400 << '\n';
  }
}

void runPercolator(std::vector<std::string> args) {
  args.insert(args.begin(), {"percolator", "-v", "0"});
  std::vector<char*> argv;
  for (std::string& arg : args) argv.push_back(&arg[0]);
  argv.push_back(nullptr);
  Caller caller;
  ASSERT_TRUE(caller.parseOptions(static_cast<int>(args.size()), argv.data()));
  ASSERT_TRUE(caller.run());
}

std::string readFile(const std::string& fileName) {
  std::ifstream in(fileName);
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

}  // namespace

// A run with protein inference used to leave process-wide settings behind:
// a later plain run drew other random numbers for the cross validation
// split, estimated its PEPs on a list taken as reversed, and then crashed
// on the missing protein estimator.
TEST(CallerTest, RunsGiveTheSameResultsAfterOtherRuns) {
  TmpDir tmpDir;  // removed with its content at the end
  fs::path dir = tmpDir.getPath();
  const std::string pin = (dir / "input.pin").string();
  writeInput(pin);
  auto file = [&](const std::string& name) { return (dir / name).string(); };
  auto plainRun = [&](const std::string& tag) {
    runPercolator({"--irls-pep", "-m", file(tag + "_psms.txt"), "-r",
                   file(tag + "_peptides.txt"), pin});
  };
  auto proteinRun = [&](const std::string& tag) {
    runPercolator({"-f", "auto", "-P", "decoy_", "-m",
                   file(tag + "_psms.txt"), "-r", file(tag + "_peptides.txt"),
                   "-l", file(tag + "_proteins.txt"), pin});
  };

  plainRun("plain1");
  proteinRun("protein1");
  plainRun("plain2");
  proteinRun("protein2");

  for (const char* suffix : {"_psms.txt", "_peptides.txt"}) {
    const std::string first = readFile(file(std::string("plain1") + suffix));
    EXPECT_FALSE(first.empty()) << suffix;
    EXPECT_EQ(first, readFile(file(std::string("plain2") + suffix))) << suffix;
  }
  for (const char* suffix : {"_psms.txt", "_peptides.txt", "_proteins.txt"}) {
    const std::string first = readFile(file(std::string("protein1") + suffix));
    EXPECT_FALSE(first.empty()) << suffix;
    EXPECT_EQ(first, readFile(file(std::string("protein2") + suffix)))
        << suffix;
  }
}

// -v and --no-terminate of one run used to stay set for every later run of
// the process: only main() restored them.
TEST(CallerTest, VerbosityAndNoTerminateDoNotCarryOver) {
  Globals::clean();  // earlier tests leave their verbosity behind
  const int defaultVerbosity = Globals::getInstance()->getVerbose();
  ASSERT_NE(0, defaultVerbosity);
  ASSERT_FALSE(Globals::getInstance()->getNoTerminate());
  {
    std::vector<std::string> args = {"percolator", "-v", "0", "--no-terminate",
                                     "missing.pin"};
    std::vector<char*> argv;
    for (std::string& arg : args) argv.push_back(&arg[0]);
    Caller caller;
    ASSERT_TRUE(caller.parseOptions(static_cast<int>(args.size()), argv.data()));
    EXPECT_EQ(0, Globals::getInstance()->getVerbose());
    EXPECT_TRUE(Globals::getInstance()->getNoTerminate());
  }
  {
    Caller caller;
    EXPECT_EQ(defaultVerbosity, Globals::getInstance()->getVerbose());
    EXPECT_FALSE(Globals::getInstance()->getNoTerminate());
  }
  Globals::clean();
}

namespace {

bool parseOnly(std::vector<std::string> args) {
  args.insert(args.begin(), {"percolator", "-v", "0"});
  std::vector<char*> argv;
  for (std::string& arg : args) argv.push_back(&arg[0]);
  argv.push_back(nullptr);
  Caller caller;
  return caller.parseOptions(static_cast<int>(args.size()), argv.data());
}

}  // namespace

// "-m x.pin x.pin" used to empty the input while the options were parsed.
TEST(CallerTest, RejectsOutputsThatOverwriteInputsOrEachOther) {
  TmpDir tmpDir;
  fs::path dir = tmpDir.getPath();
  const std::string pin = (dir / "input.pin").string();
  writeInput(pin);
  const std::string content = readFile(pin);
  const std::string out = (dir / "out.txt").string();
  const std::string sameOut = (dir / "." / "out.txt").string();

  EXPECT_THROW(parseOnly({"-m", pin, pin}), MyException);
  EXPECT_THROW(parseOnly({"-J", pin, "-m", out, pin}), MyException);
  EXPECT_THROW(parseOnly({"-w", pin, pin}), MyException);
  EXPECT_EQ(content, readFile(pin));

  EXPECT_THROW(parseOnly({"-m", out, "-r", sameOut, pin}), MyException);
  EXPECT_THROW(parseOnly({"-M", out, "-w", out, pin}), MyException);
  EXPECT_TRUE(parseOnly({"-m", out, "-r", (dir / "pep.txt").string(), pin}));
#ifndef _WIN32
  // devices may be shared
  EXPECT_TRUE(parseOnly({"-m", "/dev/null", "-r", "/dev/null", pin}));
#endif
}
