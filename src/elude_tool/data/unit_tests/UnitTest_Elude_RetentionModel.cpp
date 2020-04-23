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
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/* This file include test cases for the RTModel class */
#include <gtest/gtest.h>

#include "RetentionModel.h"
#include "PSMDescription.h"
#include "DataManager.h"
#include "Normalizer.h"
#include "Globals.h"

#define PATH_TO_DATA string("")
#define PATH_TO_WRITABLE string("")

class RetentionModelTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file = PATH_TO_DATA + "elude/standalone/train.txt";
     train_file_ptms = PATH_TO_DATA + "elude/standalone/train_2.txt";
     test_file = PATH_TO_DATA + "elude/standalone/test.txt";
     DataManager::LoadPeptides(train_file, true, true, psms, aa_alphabet);
     DataManager::LoadPeptides(train_file_ptms, true, true, psms_ptms, aa_alphabet_ptms);
     DataManager::LoadPeptides(test_file, false, true, test_psms, aa_alphabet_test);
     feature_table = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, psms);
     feature_table_ptms = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, psms_ptms);
     feature_table_test = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, test_psms);

     Normalizer::setType(Normalizer::UNI);
     rtmodel = new RetentionModel(Normalizer::getNormalizer());
     Globals::getInstance()->setVerbose(1);
   }

   virtual void TearDown() {
     delete rtmodel;
     if (feature_table) {
       delete[] feature_table;
     }
     if (feature_table_ptms) {
       delete[] feature_table_ptms;
     }
     if (feature_table_test) {
       delete[] feature_table_test;
     }
   }

   RetentionModel* rtmodel;
   string train_file, train_file_ptms, test_file;
   vector<PSMDescription> psms, psms_ptms, test_psms;
   set<string> aa_alphabet, aa_alphabet_ptms, aa_alphabet_test;
   double *feature_table, *feature_table_ptms, *feature_table_test;
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

  rtmodel->NormalizeFeatures(true, tmp);

  EXPECT_FLOAT_EQ(0.0, tmp[0].getRetentionFeatures()[0]);
  EXPECT_FLOAT_EQ(1.0, tmp[1].getRetentionFeatures()[0]);
  EXPECT_NEAR(0.17948718, tmp[2].getRetentionFeatures()[0], 0.001);

  EXPECT_FLOAT_EQ(1.0, tmp[0].getRetentionFeatures()[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[1].getRetentionFeatures()[no_features - 1]);
  EXPECT_FLOAT_EQ(0.0, tmp[2].getRetentionFeatures()[no_features - 1]);
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
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet, rf.amino_acids_alphabet(), false));
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet, rf.amino_acids_alphabet(), true));
  EXPECT_FALSE(rtmodel->IsSetIncluded(aa_alphabet_ptms, rf.amino_acids_alphabet(), false));
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet_ptms, rf.amino_acids_alphabet(), true));
  index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, true, psms_ptms);
  rf = rtmodel->retention_features();
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet_ptms, rf.amino_acids_alphabet(), false));
  EXPECT_TRUE(rtmodel->IsSetIncluded(aa_alphabet, rf.amino_acids_alphabet(), false));
}

TEST_F(RetentionModelTest, PredictRTTest) {
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, false, psms_ptms);
  rtmodel->TrainRetentionModel(aa_alphabet_ptms, index, true, psms_ptms);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet, false, "", psms));
  EXPECT_NEAR(40.3878, psms[22].getRetentionTime(), 0.01);
  EXPECT_NEAR(39.3343, psms[22].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(21.3787, psms[100].getRetentionTime(), 0.01);
  EXPECT_NEAR(23.2295, psms[100].getPredictedRetentionTime(), 0.01);
}

TEST_F(RetentionModelTest, PredictRTTestPtms) {
  Globals::getInstance()->setVerbose(1);
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet, false, psms);
  rtmodel->TrainRetentionModel(aa_alphabet, index, true, psms);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet_ptms, true, "", psms_ptms));
  EXPECT_NEAR(31.9043, psms_ptms[0].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(33.2027, psms_ptms[10].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(22.6466, psms_ptms[psms_ptms.size() - 1].getPredictedRetentionTime(), 0.01);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet_test, false, "", test_psms));
  EXPECT_NEAR(22.062, test_psms[9].getPredictedRetentionTime(), 0.01);
}

TEST_F(RetentionModelTest, SaveModelToFileTest) {
  string tmp = PATH_TO_WRITABLE + "tmp.txt";
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, false, psms_ptms);
  rtmodel->TrainRetentionModel(aa_alphabet_ptms, index, true, psms_ptms);
  rtmodel->SaveModelToFile(tmp);
  ifstream in(tmp.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE();
  }
  string line;
  getline(in, line);
  while (in && (line.find("Number_features") == string::npos)) {
     getline(in, line);
  }
  EXPECT_TRUE(line == "Number_features 44");
  getline(in, line);
  EXPECT_TRUE(line == "Active_groups 010");
  string name;
  double val;
  in >> name >> val;
  EXPECT_TRUE(name == "Sub");
  EXPECT_NEAR(41.288, val, 0.01);
  in >> name >> val;
  EXPECT_TRUE(name == "Div");
  EXPECT_NEAR(30.96125, val, 0.01);
  in >> name >> val;
  EXPECT_TRUE(name == "VSub");
  EXPECT_NEAR(-1.10973, val, 0.01);
  for(int i = 0; i < 43; ++i) {
    in >> val;
  }
  EXPECT_NEAR(0.0, val, 0.01);
  in >> name >> val;
  EXPECT_TRUE(name == "VDiv");
  EXPECT_NEAR(6.89779, val, 0.01);
  for(int i = 0; i < 42; ++i) {
    in >> val;
  }
  EXPECT_NEAR(3.0, val, 0.01);
  in >> name;
  in >> name >> val;
  EXPECT_TRUE(name == "Index");
  EXPECT_EQ(23, val);
  getline(in, line);
  getline(in, line);
  EXPECT_TRUE(line == "AA_alphabet 23 A C D E E[unimod:25] F G H I K L M N P Q R S S[unimod:21] T V W Y Y[unimod:21]");
  in.close();
  remove(tmp.c_str());
}

TEST_F(RetentionModelTest, LoadModelFromFileTest) {
  string tmp = PATH_TO_WRITABLE + "tmp.txt";
  map<string, double> index = rtmodel->BuildRetentionIndex(aa_alphabet_ptms, false, psms_ptms);
  rtmodel->TrainRetentionModel(aa_alphabet_ptms, index, true, psms_ptms);
  rtmodel->SaveModelToFile(tmp);
  Globals::getInstance()->setVerbose(1);
  delete rtmodel;
  rtmodel = new RetentionModel(Normalizer::getNormalizer());
  rtmodel->LoadModelFromFile(tmp);
  EXPECT_EQ(0,rtmodel->PredictRT(aa_alphabet, false, "", psms));
  EXPECT_NEAR(40.3878, psms[22].getRetentionTime(), 0.01);
  EXPECT_NEAR(39.3343, psms[22].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(21.3787, psms[100].getRetentionTime(), 0.01);
  remove(tmp.c_str());
}

