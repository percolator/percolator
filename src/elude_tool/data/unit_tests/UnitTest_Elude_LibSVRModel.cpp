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
/* This file include test cases for the LibSVRModel.cpp class */
#include <gtest/gtest.h>

#include "LibSVRModel.h"

#define PATH_TO_DATA string("")
#define PATH_TO_WRITABLE string("")

class LibSVRModelTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file = PATH_TO_DATA + "elude/standalone/train.txt";
     rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
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

   LibSVRModel model;
   RetentionFeatures rf;
   vector<PSMDescription> psms;
   set<string> aa_alphabet;
   string train_file;
   double *feature_table;
   int no_features;
};

TEST_F(LibSVRModelTest, TrainAndPredictBasicTest) {
  EXPECT_TRUE(feature_table != NULL) << "TrainAndPredictBasicTest error. The feature table is not initialized." << endl;
  model.setRBFSVRParam(0.01, 0.05, 5);
  model.TrainModel(psms, no_features);
  EXPECT_FALSE(model.IsModelNull()) << "TrainAndPredictBasicTest error. Null model." << endl; ;
  int len = psms.size();
  EXPECT_FLOAT_EQ(0.0, psms[len - 1].getPredictedRetentionTime());
  psms[len - 1].getPredictedRetentionTime() = model.PredictRT(no_features, psms[len - 1].getRetentionFeatures());
  // TO DO: double check that this is correct
  EXPECT_NEAR(35.5, psms[len - 1].getPredictedRetentionTime(), 0.5);
}

TEST_F(LibSVRModelTest, EstimatePredictionErrorTest) {
  vector<PSMDescription> test_psms;
  int len = psms.size();
  test_psms.push_back(psms[0]);
  test_psms.push_back(psms[len - 1]);

  model.setRBFSVRParam(0.01, 0.05, 5);
  model.TrainModel(psms, no_features);
  double pred1 =  psms[0].getRetentionTime() - model.PredictRT(no_features, psms[0].getRetentionFeatures());
  double pred2 =  psms[len - 1].getRetentionTime() - model.PredictRT(no_features, psms[len - 1].getRetentionFeatures());
  double error = model.EstimatePredictionError(no_features, test_psms);
  EXPECT_NEAR((pred1*pred1 + pred2*pred2) / 2.0, error, 0.01) << "EstimatePredictionErrorTest does not give the correct results" << endl;
}

TEST_F(LibSVRModelTest, ComputeKFoldValidationTest) {
   model.setRBFSVRParam(0.01, 0.05, 5);
  double err = model.ComputeKFoldValidation(psms, no_features);
  //TO DO: double check that this is correct
  EXPECT_NEAR(190.640, err, 1.0) << "ComputeKFoldValidationTest did not provide the correct result " << endl;
}

/*
// TO DO: check why this is SOOOO SLOW
// Check that the values are correct
TEST_F(LibSVRModelTest, CalibrateLinearModelTest) {
  model.SetSVRType(LibSVRModel::LINEAR_SVR);
  model.CalibrateLinearModel(psms, no_features);
  svm_parameter parameters = model.svr_parameters();
  EXPECT_FLOAT_EQ(0.5,parameters.C);
  EXPECT_FLOAT_EQ(0.001,parameters.p);
} */

TEST_F(LibSVRModelTest, CalibrateRBFModelTest) {
  model.CalibrateRBFModel(psms, no_features);
  svm_parameter parameters = model.svr_parameters();
  EXPECT_FLOAT_EQ(16.0,parameters.C);
  EXPECT_FLOAT_EQ(1.0,parameters.p);
  EXPECT_NEAR(0.03125,parameters.gamma, 0.001);
}


