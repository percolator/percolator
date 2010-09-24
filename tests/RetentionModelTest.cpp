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
/* This file include test cases for the RTModel class */
#include <gtest/gtest.h>

#include "RetentionModel.h"
#include "PSMDescription.h"
#include "DataManager.h"

class RetentionModelTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file = "./../bin/data/elude_test/standalone/train.txt";
     rf = rtmodel.retention_features();
     rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
     DataManager::LoadPeptides(train_file, true, true, psms, aa_alphabet);
     no_features = rf.GetTotalNumberFeatures();
     feature_table = NULL;
     feature_table = DataManager::InitFeatureTable(no_features, psms);
     rf.ComputeRetentionFeatures(psms);
   }

   virtual void TearDown() {
     if (feature_table) {
       delete[] feature_table;
     }
   }

   RetentionModel rtmodel;
   string train_file;
   RetentionFeatures rf;
   vector<PSMDescription> psms;
   set<string> aa_alphabet;
   int no_features;
   double *feature_table;
};

TEST_F(RetentionModelTest, NormalizeFeatures) {
  vector<PSMDescription> tmp;
  int last = psms.size()-1;

  tmp.push_back(psms[0]);
  tmp.push_back(psms[35]);
  tmp.push_back(psms[last]);

  rtmodel.NormalizeFeatures(tmp);

  EXPECT_FLOAT_EQ(0.0, tmp[0].retentionFeatures[0]);
  EXPECT_FLOAT_EQ(1.0, tmp[1].retentionFeatures[0]);
  EXPECT_FLOAT_EQ(0.17948718, tmp[2].retentionFeatures[0]);

  EXPECT_FLOAT_EQ(1.0, tmp[0].retentionFeatures[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[1].retentionFeatures[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[2].retentionFeatures[no_features - 1]);
}

//(std::vector<PSMDescription> &psms) {
