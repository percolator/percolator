/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@cbr.su.se>

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
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/* This file include test cases for the EludeCaller class */
#include <gtest/gtest.h>
#include <unistd.h> // for getcwd
#include <stdlib.h>// for MAX_PATH
#include "EludeCaller.h"

class EludeCallerTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file1 = "./../bin/data/elude_test/standalone/train.txt";
     test_file1 =  "./../bin/data/elude_test/standalone/test.txt";
     train_file2 = "./../bin/data/elude_test/standalone/train_1.txt";
     test_file2 = "./../bin/data/elude_test/standalone/test_1.txt";
     tmp = "./../bin/data/elude_test/standalone/tmp.txt";
     psms_.push_back(PSMDescription(10, 1));
     psms_.push_back(PSMDescription(10, 3));
     psms_.push_back(PSMDescription(10, 12));
     psms_.push_back(PSMDescription(10, 15));
     psms_.push_back(PSMDescription(8, 10));
     psms_.push_back(PSMDescription(6, 7));
     psms_.push_back(PSMDescription(10, 30));
     psms_.push_back(PSMDescription(10, 8));
     psms_.push_back(PSMDescription(10, 17));
     psms_.push_back(PSMDescription(10, 20));
     psms_.push_back(PSMDescription(10, 21));
   }

   virtual void TearDown() {
   }

   EludeCaller caller;
   string train_file1, train_file2;
   string test_file1, test_file2;
   string tmp;
   vector<PSMDescription> psms_;
};

TEST_F(EludeCallerTest, TestProcessTrainDataContext) {
  caller.set_train_file(train_file1);
  caller.set_test_file(test_file1);
  caller.set_context_format(true);
  caller.set_in_source_file(tmp);

  // no special argument
  caller.ProcessTrainData();
  vector<PSMDescription> train = caller.train_psms();
  vector<PSMDescription> test = caller.test_psms();
  EXPECT_EQ(99, train.size());
  EXPECT_EQ(1252, test.size());

  ifstream in(tmp.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE() <<  "TestProcessTrainData error: unable to open " << tmp << endl ;
  } else {
    char name[256];
    string peptide, set;
    double rt;
    in.getline(name, 256);
    in.getline(name, 256);
    in.getline(name, 256);
    in >> peptide >> rt >> set;
    EXPECT_EQ("K.REELQNVIIAQR.K", peptide);
    EXPECT_EQ(23.248, rt);
    EXPECT_EQ("train", set);
    remove(tmp.c_str());
  }

  // remove duplicates from test
  caller.set_remove_duplicates(true);
  (caller.train_psms()).clear();
  (caller.test_psms()).clear();
  caller.ProcessTrainData();
  train = caller.train_psms();
  test = caller.test_psms();
  EXPECT_EQ(99, train.size());
  EXPECT_EQ(1251, test.size());
  remove(tmp.c_str());

  (caller.train_psms()).clear();
  (caller.test_psms()).clear();
  caller.set_remove_in_source(true);
  caller.ProcessTrainData();
  EXPECT_EQ(99, caller.train_psms().size());
  EXPECT_EQ(1251, caller.test_psms().size());
  remove(tmp.c_str());

  (caller.train_psms()).clear();
   (caller.test_psms()).clear();
   caller.set_remove_in_source(true);
   caller.set_non_enzymatic(true);
   caller.ProcessTrainData();
   EXPECT_EQ(92, caller.train_psms().size());
   EXPECT_EQ(1188, caller.test_psms().size());
   remove(tmp.c_str());
}

TEST_F(EludeCallerTest, TestProcessTrainDataNoContext) {
  caller.set_train_file(train_file2);
  caller.set_test_file(test_file2);
  caller.set_in_source_file(tmp);
  caller.set_remove_common_peptides(true);
  caller.set_remove_in_source(true);
  caller.set_remove_duplicates(true);
  caller.set_test_includes_rt(true);
  caller.ProcessTrainData();
  EXPECT_EQ(135, caller.train_psms().size());
  EXPECT_EQ(52, caller.test_psms().size());
  vector<PSMDescription> psms = caller.train_psms();
  vector<PSMDescription>::iterator it = psms.begin();
  int count = 0;
  for( ; it != psms.end(); ++it)
  {
    if (it->peptide == "SNYNFEKPFLWLAR") {
      ++count;
    }
    EXPECT_FALSE("DEGWMAEHMLIMGVTRPCGR" == it->peptide);
  }
  EXPECT_EQ(1, count);
  remove(tmp.c_str());
}

TEST_F(EludeCallerTest, TestProcessTestDataNoTrain) {
  Globals::getInstance()->setVerbose(1);
  caller.set_test_file(test_file1);
  caller.set_in_source_file(tmp);
  caller.set_remove_common_peptides(true);
  caller.set_remove_in_source(true);
  caller.set_remove_duplicates(true);
  caller.set_non_enzymatic(true);
  caller.set_test_includes_rt(false);
  caller.set_context_format(true);
  caller.ProcessTestData();
  EXPECT_EQ(1188, caller.test_psms().size());
}

TEST_F(EludeCallerTest, TestTrainTestModel) {
  caller.set_train_file(train_file1);
  caller.set_test_file(test_file1);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_test_includes_rt(false);
  caller.set_context_format(true);

  Globals::getInstance()->setVerbose(1);
  caller.Run();
  EXPECT_EQ(99, caller.train_psms().size());
  vector<PSMDescription> test_psms = caller.test_psms();
  EXPECT_EQ(1252, test_psms.size());
  EXPECT_NEAR(22.496, test_psms[9].predictedTime, 0.01);
}

TEST_F(EludeCallerTest, TestComputeWindow) {
  EXPECT_NEAR(20.0, EludeCaller::ComputeWindow(psms_), 0.01);
}

TEST_F(EludeCallerTest, TestComputeRankCorrelation) {
  EXPECT_NEAR(0.4818, EludeCaller::ComputeRankCorrelation(psms_), 0.01);
}

TEST_F(EludeCallerTest, TestComputePearsonCorrelation) {
  EXPECT_NEAR(0.275, EludeCaller::ComputePearsonCorrelation(psms_), 0.01);
}
