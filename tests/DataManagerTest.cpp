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
/* This file include test cases for the DataManager class */
#include <gtest/gtest.h>

#include "DataManager.h"
#include "PSMDescription.h"
#include "Globals.h"
#include "Enzyme.h"

class DataManagerTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     train_file1 = "./../bin/data/elude_test/standalone/train.txt";
     train_file2 = "./../bin/data/elude_test/standalone/train_1.txt";
     test_file1 = "./../bin/data/elude_test/standalone/test_2.txt";
     string tmp[] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};
     basic_alphabet.insert(tmp, tmp + 20);
     Globals::getInstance()->setVerbose(5);
     Enzyme::setEnzyme(Enzyme::TRYPSIN);
   }
   virtual void TearDown() { }

   DataManager dm;
   string train_file1, train_file2, test_file1;
   set<string> basic_alphabet;
};

TEST_F(DataManagerTest, TestLoadPeptidesRTContext) {
  vector<PSMDescription> psms;
  set<string> aa_alphabet;
  // case 1: includes rt, context, but no ptms; the alphabet should be just the 20 aa
   DataManager::LoadPeptides( train_file1, true, true, psms, aa_alphabet);
  // check that the number of peptides is correct and test some of them
  EXPECT_EQ(101, psms.size()) << "TestLoadPeptidesRTContext does not give the correct results for " << train_file1 << endl;
  EXPECT_EQ("K.IIGPDADFFGELVVDAAEAVR.V", psms[32].peptide) << "TestLoadPeptidesRTContext does not give the correct results for " << train_file1 << endl ;
  EXPECT_FLOAT_EQ(62.97, psms[32].retentionTime) << "TestLoadPeptidesRTContext does not give the correct results for " << train_file1 << endl;
  EXPECT_EQ("K.QIEQGEAELEAAHTVAR.I", psms[100].peptide) << "TestLoadPeptidesRTContext does not give the correct results for " << train_file1 << endl;
  EXPECT_FLOAT_EQ(21.3787, psms[100].retentionTime) << "TestLoadPeptidesRTContext does not give the correct results for " << train_file1 << endl;
  // check the alphabet
  EXPECT_EQ(aa_alphabet.size(), basic_alphabet.size());
  set<string>::iterator it = aa_alphabet.begin();
  for( ; it != aa_alphabet.end(); ++it) {
    if (basic_alphabet.find(*it) == basic_alphabet.end())
      ADD_FAILURE() << "TestLoadPeptidesRTContext does not give the correct results for " << (*it) << endl;
  }
}

TEST_F(DataManagerTest, TestLoadPeptidesRTNoContext) {
  vector<PSMDescription> psms;
  set<string> aa_alphabet;
  // case 2: includes rt, no context, no ptms; the alphabet should be just the 20 aa
   DataManager::LoadPeptides(train_file2, true, false, psms, aa_alphabet);
  // check that the number of peptides is correct and test some of them
  EXPECT_EQ(139, psms.size()) << "TestLoadPeptidesRTNoContext does not give the correct results for " << train_file2 << endl;
  EXPECT_EQ("LTNPTYGDLNHLVSLTMSGVTTCLR", psms[32].peptide) << "TestLoadPeptidesRTNoContext does not give the correct results for " << train_file2 << endl ;
  EXPECT_FLOAT_EQ(64.7802, psms[32].retentionTime) << "TestLoadPeptidesRTNoContext does not give the correct results for " << train_file2 << endl;
  EXPECT_EQ("EIGGIFTPASVTSEEEVR", psms[138].peptide) << "TestLoadPeptidesRTNoContext does not give the correct results for " << train_file2 << endl;
  EXPECT_FLOAT_EQ(44.4893, psms[138].retentionTime) << "TestLoadPeptidesRTNoContext does not give the correct results for " << train_file2 << endl;
  // check the alphabet
  EXPECT_EQ(aa_alphabet.size(), basic_alphabet.size());
  set<string>::iterator it = aa_alphabet.begin();
  for( ; it != aa_alphabet.end(); ++it) {
    if (basic_alphabet.find(*it) == basic_alphabet.end())
      ADD_FAILURE() << "TestLoadPeptidesRTNoContext does not give the correct results for " << (*it) << endl;
  }
}

TEST_F(DataManagerTest, TestLoadPeptidesNoRTContext) {
  vector<PSMDescription> psms;
  set<string> aa_alphabet;
  // case 3: no rt, with context, with ptms; the alphabet should be the 20 aa + [PHOS]
   DataManager::LoadPeptides(test_file1, false, true, psms, aa_alphabet);
  // check that the number of peptides is correct and test some of them
  EXPECT_EQ(1251, psms.size()) << "TestLoadPeptidesNoRTContext does not give the correct results for " << test_file1 << endl;
  EXPECT_EQ("K.TMEGDCEVAYTIVQEGEK.T", psms[1250].peptide) << "TestLoadPeptidesNoRTContext does not give the correct results for " << test_file1 << endl ;
  EXPECT_FLOAT_EQ(-1.0, psms[1250].retentionTime) << "TestLoadPeptidesNoRTContext does not give the correct results for " << test_file1 << endl ;
  // check the alphabet
  basic_alphabet.insert("[PHOS]S");
  basic_alphabet.insert("[PHOS]Y");
  EXPECT_EQ(aa_alphabet.size(), basic_alphabet.size());
  set<string>::iterator it = aa_alphabet.begin();
  for( ; it != aa_alphabet.end(); ++it) {
    if (basic_alphabet.find(*it) == basic_alphabet.end())
      ADD_FAILURE() <<  "TestLoadPeptidesNoRTContext does not give the correct results for " << test_file1 << endl ;
  }
}


TEST_F(DataManagerTest, TestInitCleanFeatureTable) {
  vector<PSMDescription> psms;
  set<string> aa_alphabet;
  DataManager::LoadPeptides(test_file1, false, true, psms, aa_alphabet);
  double *feat = NULL;
  feat = dm.InitFeatureTable(10, psms);
  ASSERT_TRUE(feat != NULL) << "TestInitCleanFeatureTable (Init step) error" << endl;;
  vector<PSMDescription>::iterator it = psms.begin();
  for( ; it != psms.end(); ++it) {
    EXPECT_TRUE(it->retentionFeatures != NULL) << "TestInitCleanFeatureTable (Init step) error" << endl;
  }
  dm.CleanUpTable(psms, feat);
  for(it = psms.begin(); it != psms.end(); ++it) {
    EXPECT_TRUE(it->retentionFeatures == NULL) << "TestInitCleanFeatureTable (Clean up step) error" << endl;
  }
}

TEST_F(DataManagerTest, TestRemoveDuplicates) {
  vector<PSMDescription> psms;

  PSMDescription psm1("IAMAPEPTIDE", 10.0);
  PSMDescription psm2("PEPTIDE", 20.0);
  PSMDescription psm3("IAMAPEPTIDE", 9.0);
  PSMDescription psm4("IAMAPEPTIDE", 12.0);
  psms.push_back(psm1);
  psms.push_back(psm2);
  psms.push_back(psm3);
  psms.push_back(psm4);

  DataManager::RemoveDuplicates(psms);

  EXPECT_EQ(2, psms.size()) << "TestRemoveDuplicates error (incorrect size). " << endl;
  EXPECT_EQ(string("IAMAPEPTIDE"), psms[0].peptide) << "TestRemoveDuplicates error. " << endl;
  EXPECT_EQ(9.0, psms[0].retentionTime) << "TestRemoveDuplicates error (incorrect rt) " << endl ;
  EXPECT_EQ(string("PEPTIDE"), psms[1].peptide) << "TestRemoveDuplicates error." << endl;
}

TEST_F(DataManagerTest, TestRemoveCommonPeptides) {
  vector<PSMDescription> psms1;
  PSMDescription psm1("IAMAPEPTIDE", 10.0);
  PSMDescription psm2("PEPTIDE", 20.0);
  PSMDescription psm3("IAMAPEPTIDE", 9.0);
  psms1.push_back(psm1);
  psms1.push_back(psm2);
  vector<PSMDescription> psms2;
  PSMDescription psm4("IAMAPEPTIDE", 10.0);
  psms2.push_back(psm4);

  DataManager::RemoveCommonPeptides(psms2, psms1);

  EXPECT_EQ(1, psms1.size()) << "TestRemoveCommonPeptides error (incorrect size)." << endl;
  EXPECT_EQ(string("PEPTIDE"), psms1[0].peptide) << "TestRemoveCommonPeptides error" << endl;
  EXPECT_EQ(20.0, psms1[0].retentionTime) << "TestRemoveCommonPeptides error (incorrect rt)" << endl;
}
TEST_F(DataManagerTest, TestIsFragmentOf) {
  map<string, double> idx;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  idx["R"] = 5.0;

  // child is not included in parent
  PSMDescription child("A.AAAAYS.A", 10.0);
  PSMDescription parent("Y.YASSS.A", 20.0);
  EXPECT_FALSE(DataManager::IsFragmentOf(child, parent, 1.0, idx))
    << "TestIsFragmentOf error (child not included in parent)" << endl;

  // child included in parent and nontryptic
  child.peptide = "R.AAA.A";
  parent.peptide = "R.AAAR.A";
  EXPECT_TRUE(DataManager::IsFragmentOf(child, parent, 1.0, idx))
     << "TestIsFragmentOf error (child nontryptic" << endl;

  // child in parent, but too small difference in retention
  child.peptide = "R.AAR.A";
  EXPECT_FALSE(DataManager::IsFragmentOf(child, parent, 30.0, idx))
    << "TestIsFragmentOf error (child included in parent, small difference)" << endl;

  // child in parent, sufficient difference in retention
  EXPECT_TRUE(DataManager::IsFragmentOf(child, parent, 0.05, idx))
    << "TestIsFragmentOf error (child included in parent, large difference)" << endl;
}
