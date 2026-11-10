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
 * Unit tests for the feature normalizers.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Globals.h"
#include "MyException.h"
#include "Normalizer.h"

namespace {

// rows of two features: a normal one and one with the given values
std::vector<std::vector<double> > makeRows(const std::vector<double>& second) {
  std::vector<std::vector<double> > rows;
  for (size_t i = 0; i < second.size(); ++i) {
    rows.push_back({static_cast<double>(i), second[i]});
  }
  return rows;
}

void setSet(Normalizer& norm, std::vector<std::vector<double> >& rows) {
  std::vector<double*> features, rtFeatures;
  for (auto& row : rows) features.push_back(row.data());
  norm.setSet(features, rtFeatures, 2, 0);
}

}  // namespace

TEST(NormalizerTest, StdvNormalizerOfOrdinaryValues) {
  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
  auto rows = makeRows({1.0, 3.0, 1.0, 3.0});
  setSet(*norm, rows);
  EXPECT_DOUBLE_EQ(2.0, norm->GetVSub()[1]);
  EXPECT_DOUBLE_EQ(1.0, norm->GetVDiv()[1]);
}

// Ten copies of 0.1 add up to 0.9999999999999999, so their mean is not 0.1.
// The deviations from it gave a standard deviation of 1.4e-17, and the
// column became a constant -1 or 1, a second bias, instead of zeros.
TEST(NormalizerTest, StdvNormalizerOfAConstantColumn) {
  for (double value : {0.1, 0.7, -1e-5, 1.0 / 3.0, 3.0, 0.0}) {
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    auto rows = makeRows(std::vector<double>(10, value));
    setSet(*norm, rows);
    EXPECT_EQ(value, norm->GetVSub()[1]) << value;
    EXPECT_EQ(1.0, norm->GetVDiv()[1]) << value;
    EXPECT_EQ(0.0, norm->normalize(value, 1)) << value;
    // the other column is not affected
    EXPECT_DOUBLE_EQ(4.5, norm->GetVSub()[0]);
    EXPECT_DOUBLE_EQ(std::sqrt(8.25), norm->GetVDiv()[0]);
  }
}

// the same for a constant retention feature
TEST(NormalizerTest, StdvNormalizerOfAConstantRetentionColumn) {
  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
  auto rows = makeRows(std::vector<double>(10, 2.0));
  std::vector<std::vector<double> > rtRows;
  for (int i = 0; i < 10; ++i) rtRows.push_back({0.1, i % 2 ? 1.0 : 3.0});
  std::vector<double*> features, rtFeatures;
  for (auto& row : rows) features.push_back(row.data());
  for (auto& row : rtRows) rtFeatures.push_back(row.data());
  norm->setSet(features, rtFeatures, 2, 2);
  EXPECT_EQ(0.1, norm->GetVSub()[2]);
  EXPECT_EQ(1.0, norm->GetVDiv()[2]);
  EXPECT_DOUBLE_EQ(2.0, norm->GetVSub()[3]);
  EXPECT_DOUBLE_EQ(1.0, norm->GetVDiv()[3]);
}

// two values of 1e308 make the sum, and so the mean, infinite
TEST(NormalizerTest, StdvNormalizerRejectsAnInfiniteMean) {
  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
  auto rows = makeRows({1e308, 1e308, 1.0, 2.0});
  EXPECT_THROW(setSet(*norm, rows), MyException);
}

// one value of 1e200 has a finite mean, but its squared deviation overflows
TEST(NormalizerTest, StdvNormalizerRejectsAnInfiniteStandardDeviation) {
  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
  auto rows = makeRows({1e200, 1.0, 2.0, 3.0});
  EXPECT_THROW(setSet(*norm, rows), MyException);
}

TEST(NormalizerTest, UniNormalizerRejectsAnInfiniteRange) {
  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::UNI);
  auto rows = makeRows({1e308, -1e308, 1.0});
  EXPECT_THROW(setSet(*norm, rows), MyException);

  std::unique_ptr<Normalizer> ok = Normalizer::create(Normalizer::UNI);
  auto okRows = makeRows({1e300, -1e300, 1.0});
  EXPECT_NO_THROW(setSet(*ok, okRows));
  EXPECT_DOUBLE_EQ(2e300, ok->GetVDiv()[1]);
}

// At -v 3 the normalization factors are printed with two digits. The
// precision of cerr must be the same afterwards: the Stage lines of -v 3
// are printed with it.
TEST(NormalizerTest, StdvNormalizerRestoresThePrecisionOfCerr) {
  const int origVerbose = Globals::getInstance()->getVerbose();
  Globals::getInstance()->setVerbose(3);
  const std::streamsize origPrecision = std::cerr.precision(6);
  ::testing::internal::CaptureStderr();

  std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
  auto rows = makeRows({1.0, 3.0, 1.0, 3.0});
  setSet(*norm, rows);
  EXPECT_EQ(6, std::cerr.precision());

  std::vector<double*> features;
  for (auto& row : rows) features.push_back(row.data());
  norm->updateSet(features, 0, 2);
  EXPECT_EQ(6, std::cerr.precision());

  const std::string printed = ::testing::internal::GetCapturedStderr();
  EXPECT_NE(std::string::npos, printed.find("Normalization factors"));
  std::cerr.precision(origPrecision);
  Globals::getInstance()->setVerbose(origVerbose);
}
