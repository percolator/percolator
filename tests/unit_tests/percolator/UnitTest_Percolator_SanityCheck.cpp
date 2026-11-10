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
 * Unit tests for reading a weights file (--init-weights).
 */

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

#include "DataSet.h"
#include "MyException.h"
#include "SanityCheck.h"

class ReadWeightsTest : public ::testing::Test {
 protected:
  void SetUp() override {
    DataSet::resetFeatureNames();
    DataSet::getFeatureNames().insertFeature("xcorr");
    DataSet::getFeatureNames().insertFeature("deltCn");
    FeatureNames::setNumFeatures(2);
    check.setInitWeightFN("w.txt");
  }
  void TearDown() override { DataSet::resetFeatureNames(); }
  std::vector<double> read(const std::string& text) {
    std::istringstream in(text);
    std::vector<double> w(3, 0.0);
    check.readWeights(in, w);
    return w;
  }
  SanityCheck check;
};

// the format written by -w: names, normalized and raw weights per split
TEST_F(ReadWeightsTest, ReadsTheRawWeightsOfTheFirstSplit) {
  EXPECT_EQ((std::vector<double>{1.5, -2.0, 0.25}),
            read("# comment\n"
                 "xcorr\tdeltCn\tm0\n"
                 "0.1\t0.2\t0.3\n"
                 "1.5\t-2\t0.25\n"
                 "xcorr\tdeltCn\tm0\n"
                 "9\t9\t9\n"
                 "9\t9\t9\n"));
  // without names, only the number of weights can be checked
  EXPECT_EQ((std::vector<double>{1.0, 2.0, 3.0}),
            read("0.1\t0.2\t0.3\n1\t2\t3\r\n"));
}

TEST_F(ReadWeightsTest, RejectsWeightsOfOtherFeatures) {
  // other order, fewer and more features
  EXPECT_THROW(read("deltCn\txcorr\tm0\n0.1\t0.2\t0.3\n1\t2\t3\n"),
               MyException);
  EXPECT_THROW(read("xcorr\tm0\n0.1\t0.2\n1\t2\n"), MyException);
  EXPECT_THROW(read("xcorr\tdeltCn\tsp\tm0\n0\t0\t0\t0\n1\t2\t3\t4\n"),
               MyException);
  EXPECT_THROW(read("0.1\t0.2\n1\t2\n"), MyException);
  EXPECT_THROW(read("0.1\t0.2\t0.3\t0.4\n1\t2\t3\t4\n"), MyException);
}

TEST_F(ReadWeightsTest, RejectsMissingOrNonFiniteWeights) {
  EXPECT_THROW(read("xcorr\tdeltCn\tm0\n0.1\t0.2\t0.3\n"), MyException);
  EXPECT_THROW(read("# only comments\n"), MyException);
  EXPECT_THROW(read("xcorr\tdeltCn\tm0\n0.1\t0.2\t0.3\n1\tnan\t3\n"),
               MyException);
  EXPECT_THROW(read("xcorr\tdeltCn\tm0\n0.1\t0.2\t0.3\n1\tx\t3\n"),
               MyException);
}

// Feature weights that are all 0, whatever m0 is, score every PSM alike;
// before, such a run ended where the bins could not be rescaled. The error
// names the file. A single feature weight that is not 0 is enough.
TEST_F(ReadWeightsTest, RejectsWeightsThatAreAllZero) {
  for (const char* raw : {"0\t0\t0\n", "0\t-0\t0.25\n", "-0\t0\t-3\n"}) {
    SCOPED_TRACE(raw);
    try {
      read(std::string("xcorr\tdeltCn\tm0\n0.1\t0.2\t0.3\n") + raw);
      ADD_FAILURE() << "no exception";
    } catch (const MyException& e) {
      const std::string message = e.what();
      EXPECT_NE(std::string::npos, message.find("weights file w.txt"));
      EXPECT_NE(std::string::npos, message.find("weight 0"));
    }
  }
  EXPECT_EQ((std::vector<double>{0.0, 1e-300, 0.0}),
            read("xcorr\tdeltCn\tm0\n0\t0\t0\n0\t1e-300\t0\n"));
}
