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
#include "Normalizer.h"

class RetentionModelTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file = "./../bin/data/elude_test/standalone/train.txt";
     train_file_ptms = "./../bin/data/elude_test/standalone/train_2.txt";
     DataManager::LoadPeptides(train_file, true, true, psms, aa_alphabet);
     DataManager::LoadPeptides(train_file_ptms, true, true, psms_ptms, aa_alphabet_ptms);
     feature_table = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, psms);
     feature_table_ptms = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, psms_ptms);
     rtmodel = new RetentionModel(Normalizer::getNormalizer());
   }

   virtual void TearDown() {
     delete rtmodel;
     if (feature_table) {
       delete[] feature_table;
     }
     if (feature_table_ptms) {
       delete[] feature_table_ptms;
     }
   }

   RetentionModel* rtmodel;
   string train_file, train_file_ptms;
   vector<PSMDescription> psms, psms_ptms;
   set<string> aa_alphabet, aa_alphabet_ptms;
   double *feature_table, *feature_table_ptms;
 };

TEST_F(RetentionModelTest, NormalizeFeaturesTest) {
  RetentionFeatures rf;
  int no_features;
  rf = rtmodel->retention_features();
  rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
  no_features = rf.GetTotalNumberFeatures();
  rf.ComputeRetentionFeatures(psms);

  vector<PSMDescription> tmp;
  int last = psms.size()-1;

  tmp.push_back(psms[0]);
  tmp.push_back(psms[35]);
  tmp.push_back(psms[last]);

  rtmodel->NormalizeFeatures(tmp);

  EXPECT_FLOAT_EQ(0.0, tmp[0].retentionFeatures[0]);
  EXPECT_FLOAT_EQ(1.0, tmp[1].retentionFeatures[0]);
  EXPECT_NEAR(0.17948718, tmp[2].retentionFeatures[0], 0.001);

  EXPECT_FLOAT_EQ(1.0, tmp[0].retentionFeatures[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[1].retentionFeatures[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[2].retentionFeatures[no_features - 1]);
}

TEST_F(RetentionModelTest, BuildRetentionIndexNoPtmsTest) {
  PSMDescription::setPSMSet(psms);
  PSMDescription::normalizeRetentionTimes(psms);
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet, true, psms);
  EXPECT_NEAR(0.0194839, index["A"], 0.01);
  EXPECT_NEAR(0.884652, index["L"], 0.01);
  EXPECT_NEAR(0.193366, index["Y"], 0.01);
}

TEST_F(RetentionModelTest, BuildRetentionIndexTestNoPtmsTest) {
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, false, psms_ptms);
  EXPECT_NEAR(0.075016, index["E[unimod:25]"], 0.1);
  EXPECT_NEAR(0.0295835, index["A"], 0.1);
  EXPECT_NEAR(0.112165, index["S[unimod:21]"], 0.1);
}

TEST_F(RetentionModelTest, TrainRetentionModelNoPtmsTest) {
  EXPECT_TRUE(rtmodel->IsModelNull());
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet, false, psms);
  rtmodel->TrainRetentionModel(aa_alphabet, index, true, psms);
  EXPECT_FALSE(rtmodel->IsModelNull());
}

TEST_F(RetentionModelTest, TrainRetentionModelPtmsTest) {
  EXPECT_TRUE(rtmodel->IsModelNull());
  PSMDescription::setPSMSet(psms_ptms);
  PSMDescription::normalizeRetentionTimes(psms_ptms);
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, true, psms_ptms);
  rtmodel->TrainRetentionModel(aa_alphabet_ptms, index, true, psms_ptms);
  EXPECT_FALSE(rtmodel->IsModelNull());
}

TEST_F(RetentionModelTest, IsSetIncludedTest) {
  PSMDescription::setPSMSet(psms);
  PSMDescription::normalizeRetentionTimes(psms);
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet, true, psms);
  RetentionFeatures rf = rtmodel->retention_features();
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet, rf.amino_acids_alphabet()));
  EXPECT_FALSE(rtmodel->IsSetIncluded(aa_alphabet_ptms, rf.amino_acids_alphabet()));
  index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, true, psms_ptms);
  rf = rtmodel->retention_features();
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet_ptms, rf.amino_acids_alphabet()));
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet, rf.amino_acids_alphabet()));
}

TEST_F(RetentionModelTest, PredictRTTest) {
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, false, psms_ptms);
  rtmodel->TrainRetentionModel(aa_alphabet_ptms, index, true, psms_ptms);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet, psms));
  EXPECT_NEAR(40.3878, psms[22].retentionTime, 0.01);
  EXPECT_NEAR(39.3343, psms[22].predictedTime, 0.01);
  EXPECT_NEAR(21.3787, psms[100].retentionTime, 0.01);
  EXPECT_NEAR(23.2295, psms[100].predictedTime, 0.01);
}

/*
TEST_F(RetentionModelTest, PredictRTTestPtms) {
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet, false, psms);
  rtmodel->TrainRetentionModel(aa_alphabet, index, true, psms);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet_ptms, psms_ptms));
}*/
