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
/* This file include test cases for LibsvmWrapper */
#include <gtest/gtest.h>

#include "LibsvmWrapper.h"
#include "PSMDescription.h"
#include "DataManager.h"
#include "RetentionFeatures.h"

class LibsvmWrapperTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file = "./../bin/data/elude_test/standalone/train.txt";
     DataManager::LoadPeptides(train_file, true, true, psms, aa_alphabet);
     no_features = rf.GetTotalNumberFeatures();
     feature_table = DataManager::InitFeatureTable(no_features, psms);
     rf.ComputeRetentionFeatures(psms);
     Globals::getInstance()->setVerbose(1);
   }

   virtual void TearDown() {
     DataManager::CleanUpTable(psms, feature_table);
     feature_table = NULL;
   }

   RetentionFeatures rf;
   vector<PSMDescription> psms;
   set<string> aa_alphabet;
   string train_file;
   double *feature_table;
   int no_features;
};


TEST_F(LibsvmWrapperTest, SetUp) {

  EXPECT_EQ(0.0, 0.0);
}

