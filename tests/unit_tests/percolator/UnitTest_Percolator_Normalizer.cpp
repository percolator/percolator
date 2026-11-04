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

#include <memory>
#include <vector>

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
