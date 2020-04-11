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
/* This file include test cases for the EludeCaller class */
#include <gtest/gtest.h>
#include <unistd.h> // for getcwd
#include <stdlib.h>// for MAX_PATH
#include <algorithm>
#include <fstream>
#include "EludeCaller.h"
#include "Globals.h"

#define PATH_TO_DATA string("")
#define PATH_TO_WRITABLE string("")

class EludeCallerTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file1 = PATH_TO_DATA + "elude/standalone/train.txt";
     test_file1 =  PATH_TO_DATA + "elude/standalone/test.txt";
     train_file2 = PATH_TO_DATA + "elude/standalone/train_1.txt";
     test_file2 = PATH_TO_DATA + "elude/standalone/test_1.txt";
     tmp = PATH_TO_WRITABLE + "tmp.txt";
     calibration_file = PATH_TO_DATA + "elude/calibrate_data/calibrate.txt";
     lib_path = PATH_TO_DATA + "elude/calibrate_data/test_lib";
     test_calibration =  PATH_TO_DATA + "elude/calibrate_data/test.txt";
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
     Globals::getInstance()->setVerbose(1);
   }

   virtual void TearDown() {
   }
   EludeCaller caller;
   string train_file1, train_file2;
   string test_file1, test_file2;
   string calibration_file, lib_path, test_calibration;
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

TEST_F(EludeCallerTest, TestProcessTrainDataNoContext)  {
  caller.clear_train_psms();
  caller.clear_test_psms();
  caller.set_train_file(train_file2);
  caller.set_test_file(test_file2);
  caller.set_in_source_file(tmp);
  caller.set_remove_common_peptides(true);
  caller.set_remove_in_source(true);
  caller.set_remove_duplicates(true);
  caller.set_test_includes_rt(true);
  caller.ProcessTrainData();
  EXPECT_EQ(135, caller.train_psms().size());
  EXPECT_EQ(53, caller.test_psms().size());
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

  caller.Run();
  EXPECT_EQ(99, caller.train_psms().size());
  vector<PSMDescription> test_psms = caller.test_psms();
  EXPECT_EQ(1252, test_psms.size());
  sort(test_psms.begin(), test_psms.end());
  EXPECT_NEAR(28.3021, test_psms[9].getPredictedRetentionTime(), 2.0);
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

TEST_F(EludeCallerTest, TestSaveLoadModel) {
  caller.set_train_file(train_file1);
  caller.set_save_model_file(tmp);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);

  caller.Run();

  EludeCaller caller2;
  caller2.set_load_model_file(tmp);
  caller2.set_test_file(test_file1);
  caller2.set_remove_common_peptides(false);
  caller2.set_remove_in_source(false);
  caller2.set_remove_duplicates(false);
  caller2.set_non_enzymatic(false);
  caller2.set_context_format(true);
  caller2.Run();
  vector<PSMDescription> test_psms = caller2.test_psms();
  EXPECT_EQ(1252, test_psms.size());
  sort(test_psms.begin(), test_psms.end());
  EXPECT_NEAR(29.6867, test_psms[9].getPredictedRetentionTime(), 2.0);
  remove(tmp.c_str());
}

TEST_F(EludeCallerTest, TestListDirFiles) {
  string dir = PATH_TO_DATA + "elude/standalone/";
  vector<string> files = EludeCaller::ListDirFiles(dir);

  EXPECT_EQ(8.0, files.size());
  sort(files.begin(), files.end());
  int pos = files[0].rfind("/");
  string f = files[0].substr(pos+1, files[0].length() - pos);
  EXPECT_TRUE("test.txt" == f);
  pos = files[files.size() - 1].rfind("/");
  f = files[files.size() - 1].substr(pos+1, files[files.size() - 1].length() - pos);
  EXPECT_TRUE("train_3.txt" == f);
}

TEST_F(EludeCallerTest, TestAutomaticModelSelection) {
  caller.set_train_file(calibration_file);
  caller.set_test_file(test_calibration);
  caller.set_automatic_model_sel(true);
  caller.set_context_format(true);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_test_includes_rt(true);
  EludeCaller::set_lib_path(lib_path);
  caller.set_linear_calibration(false);
  caller.Run();

  vector<PSMDescription> test_psms = caller.test_psms();
  EXPECT_EQ(1740, test_psms.size());
  sort(test_psms.begin(), test_psms.end());
  EXPECT_NEAR(51.6007, test_psms[0].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(27.528, test_psms[1000].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(39.1311, test_psms[1739].getPredictedRetentionTime(), 0.01);
}

TEST_F(EludeCallerTest, TestFindLeastSquaresSolution) {
  vector<PSMDescription> psms2;
  psms2.push_back(PSMDescription(3, 1));
  psms2.push_back(PSMDescription(5, 2));
  psms2.push_back(PSMDescription(7, 3));
  double a = 0.0, b = 0.0;
  EludeCaller::FindLeastSquaresSolution(psms2, a, b);
  EXPECT_NEAR(2.0, a, 0.01);
  EXPECT_NEAR(1.0, b, 0.01);
}

TEST_F(EludeCallerTest, TestAutomaticModelSelectionWithCalibration) {
  caller.set_train_file(calibration_file);
  caller.set_test_file(test_calibration);
  caller.set_automatic_model_sel(true);
  caller.set_context_format(true);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_test_includes_rt(true);
  EludeCaller::set_lib_path(lib_path);
  double cov = 1.0;
  LTSRegression::setCoverage(cov);
  caller.Run();
  pair<double, double> coeff = caller.lts_coefficients();
  vector<PSMDescription> train_psms = caller.train_psms();
  double a = 0.0, b = 0.0;
  EludeCaller::FindLeastSquaresSolution(train_psms, a, b);
  EXPECT_NEAR(a, coeff.first, 0.01);
  EXPECT_NEAR(b, coeff.second, 0.01);
  vector<PSMDescription> test_psms = caller.test_psms();
  EXPECT_EQ(1740, test_psms.size());
  sort(test_psms.begin(), test_psms.end());
  EXPECT_NEAR(51.6007 * a + b, test_psms[0].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(27.528 * a + b, test_psms[1000].getPredictedRetentionTime(), 0.01);
  EXPECT_NEAR(39.1311 * a + b , test_psms[1739].getPredictedRetentionTime(), 0.01);
}

TEST_F(EludeCallerTest, TestGetFileName) {
  string path = "D:/dir/file.txt";
  EXPECT_TRUE("file" == EludeCaller::GetFileName(path));
  path = "D:\\dir\\file";
  EXPECT_TRUE("file" == EludeCaller::GetFileName(path));
  path = "file.txt";
  EXPECT_TRUE("file" == EludeCaller::GetFileName(path));
  path = "file";
  EXPECT_TRUE("file" == EludeCaller::GetFileName(path));
  EXPECT_TRUE("train" == EludeCaller::GetFileName(train_file1));
}

TEST_F(EludeCallerTest, TestAddModelLibraryTrain) {
  caller.set_train_file(train_file1);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);
  caller.set_lib_path(PATH_TO_DATA + "elude/calibrate_data/test_lib");
  caller.set_append_model(true);
  caller.Run();

  // check that the file exists
  string file_name = PATH_TO_DATA + "elude/calibrate_data/test_lib/train.model";
  ifstream in(file_name.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE() <<  "TestAddModelLibrary error: unable to open " << file_name << endl;
  }
  in.close();
  remove(file_name.c_str());
}

TEST_F(EludeCallerTest, TestAddModelLibrarySave) {
  caller.set_train_file(train_file1);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);
  string model_file = PATH_TO_DATA + "elude/calibrate_data/test_lib/test.model";
  caller.set_save_model_file(model_file);
  caller.set_lib_path(PATH_TO_DATA + "elude/calibrate_data/test_lib");
  caller.set_append_model(true);
  caller.Run();

  // check that the file exists
  ifstream in(model_file.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE() <<  "TestAddModelLibrarySave error: unable to open " << model_file << endl;
  }
  in.close();
  remove(model_file.c_str());
}

TEST_F(EludeCallerTest, TestSaveIndexToFileTrain) {
  caller.set_train_file(train_file1);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);
  string tmp_file =  PATH_TO_WRITABLE  + "tmp";
  caller.set_index_file(tmp_file);
  caller.Run();
  ifstream in(tmp_file.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE() <<  "TestSaveIndexToFileTrain error: unable to open " <<tmp_file << endl;
  }
  string tmp;
  getline(in, tmp);
  getline(in, tmp);
  string aa;
  double val;
  for(int i = 0; i < 10; ++i) {
    in >> aa >> tmp >> val;
  }
  EXPECT_TRUE(aa == "L");
  EXPECT_NEAR(val, 0.928439, 0.5);

  in.close();
  remove(tmp_file.c_str());
}

TEST_F(EludeCallerTest, TestSaveIndexToFileAutSel) {
  caller.set_train_file(calibration_file);
  caller.set_test_file(test_calibration);
  caller.set_automatic_model_sel(true);
  caller.set_context_format(true);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_test_includes_rt(true);
  EludeCaller::set_lib_path(lib_path);
  string tmp_file = PATH_TO_WRITABLE  + "tmp";
  caller.set_index_file(tmp_file);
  double cov = 1.0;
  LTSRegression::setCoverage(cov);
  caller.Run();
  ifstream in(tmp_file.c_str(), ios::in);
  if (in.fail()) {
    ADD_FAILURE() <<  "TestSaveIndexToFileAutSel error: unable to open " <<tmp_file << endl;
  }
  string tmp;
  for(int i = 0; i < 12; ++i) {
    getline(in, tmp);
  }
  in.close();
  EXPECT_TRUE("L : 1.49573" == tmp);
  remove(tmp_file.c_str());
} 

/****************** TO BE REMOVED *******************/

/* Real experiment */
/*
TEST_F(EludeCallerTest, TestAUtomaticSelectioNSystem) {
  Globals::getInstance()->setVerbose(5);
  caller.set_train_file("/scratch/lumi_work/projects/retention_time/results/article_revised_results/data/txt_pep_0_01/clean_train/Jupiter_120_clean_train.txt");
  caller.set_test_file("/scratch/lumi_work/projects/retention_time/results/article_revised_results/data/txt_pep_0_01/clean_test/Jupiter_120_clean_test.txt");
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);
  caller.set_test_includes_rt(true);

  caller.set_automatic_model_sel(true);
  caller.set_linear_calibration(true);
  double cov = 0.95;
  LTSRegression::setCoverage(cov);
  caller.set_lib_path("/scratch/lumi_work/projects/elude_ptms/src/bin/data/elude_test/calibrate_data/test_lib");
  //caller.set_lib_path("/scratch/lumi_work/projects/elude_ptms/src/bin/data/elude_lib");
  caller.Run();
}*/

/*
TEST_F(EludeCallerTest, TestLoadModel) {
  Globals::getInstance()->setVerbose(1);
  caller.set_load_model_file("/scratch/lumi_work/projects/elude_ptms/src/bin/data/elude_test/calibrate_data/test_lib/Jupiter_60.model");
  caller.set_test_file(test_calibration);
  caller.set_remove_common_peptides(false);
  caller.set_remove_in_source(false);
  caller.set_remove_duplicates(false);
  caller.set_non_enzymatic(false);
  caller.set_context_format(true);
  caller.set_test_includes_rt(true);
  caller.Run();
  vector<PSMDescription> test_psms = caller.test_psms();
  sort(test_psms.begin(), test_psms.end());
  EXPECT_EQ(1740, test_psms.size());
  cout << test_psms[0].peptide << " " << test_psms[0].getPredictedRetentionTime() << endl;
  cout << test_psms[1000].peptide << " " << test_psms[1000].getPredictedRetentionTime() << endl;
  cout << test_psms[1739].peptide << " " << test_psms[1739].getPredictedRetentionTime() << endl;
}*/


/*

Performance measures for the test data:
  Pearson's correlation r = 0.955049
  Spearman's rank correlation rho = 0.962522
  Delta_t 95% = 23.7204 */
/*
TEST_F(EludeCallerTest, TestRunTrainTest) {
  caller.set_train_file("/scratch/lumi_work/projects/retention_time/results/article_revised_results/data/txt_pep_0_01/clean_train/Luna_240_clean_train.txt");
  caller.set_save_model_file("/scratch/lumi_work/projects/elude_ptms/src/bin/data/elude_lib/Luna_240.model");
  caller.set_test_file("/scratch/lumi_work/projects/retention_time/results/article_revised_results/data/txt_pep_0_01/clean_test/Luna_240_clean_test.txt");
  caller.set_remove_common_peptides(true);
  caller.set_remove_in_source(true);
  caller.set_remove_duplicates(true);
  caller.set_non_enzymatic(true);
  caller.set_test_includes_rt(true);
  caller.set_context_format(true);

  Globals::getInstance()->setVerbose(5);
  caller.Run();
}
*/
