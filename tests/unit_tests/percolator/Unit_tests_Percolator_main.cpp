/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@scilifelab.se>

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
 * @ Created by Mattia Tomasoni
 * Feb, 2011
 */
/*
 * Main file for Google Test
 * Just run all tests
 */
#include <gtest/gtest.h>

#include "Caller.h"

namespace {

// Every test starts from the initial process-wide settings. A run of
// Caller leaves its settings behind (protein inference after -f, the
// verbosity), and a test that failed part-way does not restore what it
// set, so without this the result of a test depended on the tests before
// it, and on their order with --gtest_shuffle.
class ResetProcessSettings : public ::testing::EmptyTestEventListener {
  void OnTestStart(const ::testing::TestInfo&) override {
    Caller::resetProcessSettings();
  }
};

}  // namespace

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  // the listeners own the object
  ::testing::UnitTest::GetInstance()->listeners().Append(
      new ResetProcessSettings);
  return RUN_ALL_TESTS();
}
