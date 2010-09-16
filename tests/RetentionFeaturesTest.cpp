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
/* This file include test cases for the RetentionFeatures class */
#include <gtest/gtest.h>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <bitset>

#include "RetentionFeatures.h"
#include "PSMDescription.h"

//using namespace std;

class RetentionFeaturesTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     PSMDescription psm1(string("PPPP[PHOS]AAAA[GLYC]"), 10.5);
     PSMDescription psm2(string("AA[PHOS]PP[GLYC]"), 11.5);
     psm1.retentionFeatures = new double[60];
     psm2.retentionFeatures = new double[60];
     psms.push_back(psm1);
     psms.push_back(psm2);
   }

   virtual void TearDown() {
     for (int i = 0; i < psms.size(); ++i)
       delete[] psms[i].retentionFeatures;
   }

   vector<PSMDescription> psms;
   RetentionFeatures rf;
};

TEST_F(RetentionFeaturesTest, TestGetIndexValue) {
  /* simple amino acids */
  double val = RetentionFeatures::GetIndexValue("S", RetentionFeatures::k_kyte_doolittle());
  EXPECT_EQ(-0.8, val) << "GetIndexValue does not give the correct results for the entry \"S\"" << endl;

  /* modified amino acid not included in index */
  val = RetentionFeatures::GetIndexValue("[PHOS]Y", RetentionFeatures::k_bulkiness());
  EXPECT_EQ(18.03, val) << "GetIndexValue does not give the correct results for the entry \"[PHOS]Y\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 6.5;
  idx["A"] = 10.0;
  val = RetentionFeatures::GetIndexValue("[PHOS]A", idx);
  EXPECT_EQ(6.5, val) << "GetIndexValue does not give the correct results for the entry \"[PHOS]A\" "
                      << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexSum) {
  /* the index does not include the modified peptide */
  string peptide = "ASA[PHOS]AY";
  double hydrophobicity_sum = RetentionFeatures::IndexSum(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(3.3, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  hydrophobicity_sum =  RetentionFeatures::IndexSum(peptide, idx);
  EXPECT_FLOAT_EQ(7.7, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\" "
                                            << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestPeptideLength) {
  /* no modified aa */
  double val = RetentionFeatures::PeptideLength("AASSY");
  EXPECT_FLOAT_EQ(5.0, val) << "PeptideLength does not calculate the correct length of the peptide \"AASSY\"" << endl;

  /* modified aa */
  val = RetentionFeatures::PeptideLength("[PHOS]AS[GLYC]Y");
  EXPECT_FLOAT_EQ(3.0, val) << "PeptideLength does not calculate the correct length of the peptide \"[PHOS]A[GLYC]Y\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexAvg) {
  /* the index does not include the modified peptide */
  string peptide = "ASA[PHOS]AY";
  double hydrophobicity_avg =  RetentionFeatures::IndexAvg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(3.3/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  peptide = "A[GLYC]SA[PHOS]AY";
  hydrophobicity_avg =  RetentionFeatures::IndexAvg(peptide, idx);
  EXPECT_FLOAT_EQ(7.7/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"A[GLYC]SA[PHOS]AY\" "
                                            << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestGetAminoAcids) {
  string peptide = "[PHOS]AAA[PHOS]Y[PHOS]Y";
  string exp_peptides[5] = {"[PHOS]A", "A", "A", "[PHOS]Y", "[PHOS]Y"};

  vector<string> res_peptides =  RetentionFeatures::GetAminoAcids(peptide);
  EXPECT_EQ(res_peptides.size(), 5);
  for(int i = 0; i < res_peptides.size(); ++i)
    EXPECT_EQ(exp_peptides[i], res_peptides[i]);
}

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbour) {
  // get the polar aa
  set<string> polar_aa = (RetentionFeatures::GetExtremeRetentionAA(RetentionFeatures::k_kyte_doolittle())).first;
  // no modified aa
  string peptide = "KAYYKYCNCYYYRA";
  double sum_neighbour_pos =  RetentionFeatures::IndexNearestNeigbour(peptide, RetentionFeatures::k_kyte_doolittle(), polar_aa);
  EXPECT_FLOAT_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbour does not give the correct results for the entry \"KAYYKYCNCYYYRA\"" << endl;

  // modified aa is in the index
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["[GLYC]C"] = -1.0;
  idx["[PHOS]Y"] = 2.0;
  idx["A"] = 5.0;
  peptide = "KAA[PHOS]AA[GLYC]CA";
  polar_aa = RetentionFeatures::GetExtremeRetentionAA(idx).first;
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbour(peptide, idx, polar_aa);
  EXPECT_FLOAT_EQ(10.0, sum_neighbour_pos) << "IndexNearestNeigbour does not give the correct results for the entry \"KAA[PHOS]AA[GLYC]CA\"" << endl;
}

/*
TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourPos) {
  // no modified aa
  string peptide = "KAYYKYCRCYYYRA";
  double sum_neighbour_pos =  RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAYYKYCRCYYYRA\"" << endl;

  // modified aa; the index does not include the modified form
  peptide = "[PHOS]K[PHOS]AYYKY[PHOS]CRCYYYR[PHOS]A";
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"[PHOS]K[PHOS]AYYKY[PHOS]CRCYYYR[PHOS]A\"" << endl;

  // modified aa is in the index
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["[GLYC]C"] = -1.0;
  idx["[PHOS]Y"] = 2.0;
  idx["A"] = 5.0;
  peptide = "KAA[PHOS]AR[GLYC]CKA";
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, idx);
  EXPECT_FLOAT_EQ(11.0, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAA[PHOS]AR[GLYC]CKA\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourNeg) {
  // no modified aa
  string peptide = "DAYYDYCECYYYEA";
  double sum_neighbour_neg =  RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAYYDYCECYYYEA\"" << endl;

  // modified aa; the index does not include the modified form
  peptide = "[PHOS]D[PHOS]AYYDY[PHOS]CECYYYE[PHOS]A";
  sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"[PHOS]D[PHOS]AYYDY[PHOS]CECYYYE[PHOS]A\"" << endl;

  // modified aa is in the index
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["[GLYC]C"] = -1.0;
  idx["[PHOS]Y"] = 2.0;
  idx["A"] = 5.0;
  peptide = "DAA[PHOS]AE[GLYC]CDA";
  sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, idx);
  EXPECT_FLOAT_EQ(11.0, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAA[PHOS]AE[GLYC]CDA\"" << endl;
}*/

TEST_F(RetentionFeaturesTest, TestIndexN) {
  /* no modified aa */
  string peptide = "DAYY";
  double hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-3.5, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"DAYY\"" << endl;

  /* second aa is modified but not included in the index */
  peptide = "Y[PHOS]DAYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-1.3, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"Y[PHOS]DAYY\"" << endl;

  /* N term is modified but not included in the index */
  peptide = "[PHOS]DAYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-3.5, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"[PHOS]DAYY\"" << endl;

  /* N term is modified and included in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "[PHOS]AE";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, idx);
  EXPECT_FLOAT_EQ(1.0, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"[PHOS]AE\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexC) {
  /* no modified aa */
  string peptide = "DAYY";
  double hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-1.3, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"DAYY\"" << endl;

  /* one but last is modified but not included in the index */
  peptide = "Y[PHOS]DAYY[PHOS]DD";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-3.5, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"Y[PHOS]DAYY[PHOS]DD\"" << endl;

  /* C term is modified but not included in the index */
  peptide = "[PHOS]DAY[PHOS]D";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(-3.5, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"[PHOS]DAY[PHOS]D\"" << endl;

  /* C term is modified and included in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "E[PHOS]A";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, idx);
  EXPECT_FLOAT_EQ(1.0, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"E[PHOS]A\"" << endl;
}

/*
TEST_F(RetentionFeaturesTest, TestIndexNC) {
  // no modified aa
  string peptide = "DAYY";
  double hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(0.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"DAYY\"" << endl;

  // second and last aa are modified but not included in the index
  peptide = "I[PHOS]DAYY[PHOS]F";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"I[PHOS]DAYY[PHOS]F\"" << endl;

  // both c and n terminus are modified but not included in the index
  peptide = "[PHOS]IAY[PHOS]F";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_FLOAT_EQ(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[PHOS]IAY[PHOS]F\"" << endl;

  // C term and N term are modified and included in the index
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "[PHOS]A";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, idx);
  EXPECT_FLOAT_EQ(1.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[PHOS]A\"" << endl;
}*/

TEST_F(RetentionFeaturesTest, TestIndexMaxPartialSum) {
  /* no modified aa */
  string peptide = "DDDFFFDDDIIDDD";
  double max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_FLOAT_EQ(9.0, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  /* modified aa but not in the index */
  peptide = "DDDA[PHOS]AADDDI[PHOS]IWDD";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_FLOAT_EQ(8.1, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;

  /* modified aa that are in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  peptide = "[PHOS]AADD[PHOS]D";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, idx, 3);
  EXPECT_FLOAT_EQ(11.0, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinPartialSum) {
  /* no modified aa */
  string peptide = "DDDFFFDDDIIDDD";
  double min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_FLOAT_EQ(-7.0, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  /* modified aa but not in the index */
  peptide = "[PHOS]D[PHOS]DDADD";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_FLOAT_EQ(-10.5, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;

  /* modified aa that are in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  peptide = "[PHOS]AADD[PHOS]D";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, idx, 3);
  EXPECT_FLOAT_EQ(8.0, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestAvgHydrophobicityIndex) {
  /* no modified aa */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  double avg_hydrophobicity = RetentionFeatures:: AvgHydrophobicityIndex(idx);
  EXPECT_FLOAT_EQ(2.75, avg_hydrophobicity) << "TestAvgHydrophobicityIndex does not give the correct result" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMaxHydrophobicSideHelix) {
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  /* peptide smaller than 9 aa */
  string peptide = "DDDDDD";
  double max_hydrophobic_side = RetentionFeatures::IndexMaxHydrophobicSideHelix(peptide, idx);
  EXPECT_FLOAT_EQ(7.841237327, max_hydrophobic_side) << "TestIndexMaxHydrophobicSideHelix does not give the correct results for the entry \"DDDDDD\"" << endl;

  /* normal peptides (with modified aa) */
  peptide = "[PHOS]ADDAA[PHOS]D[PHOS]ADDY";
  max_hydrophobic_side = RetentionFeatures::IndexMaxHydrophobicSideHelix(peptide, idx);
  EXPECT_FLOAT_EQ(11.064177772, max_hydrophobic_side) << "TestIndexMaxHydrophobicSideHelix does not give the correct results for the entry \"[PHOS]ADDAA[PHOS]D[PHOS]ADDY\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinHydrophobicSideHelix) {
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  /* peptide smaller than 9 aa */
  string peptide = "DDDDDD";
  double min_hydrophobic_side = RetentionFeatures::IndexMinHydrophobicSideHelix(peptide, idx);
  EXPECT_FLOAT_EQ(7.841237327, min_hydrophobic_side) << "TestIndexMinHydrophobicSideHelix does not give the correct results for the entry \"DDDDDD\"" << endl;

  /* normal peptides (with modified aa) */
  peptide = "[PHOS]ADDAA[PHOS]D[PHOS]ADDY";
  min_hydrophobic_side = RetentionFeatures::IndexMinHydrophobicSideHelix(peptide, idx);
  EXPECT_FLOAT_EQ(6.438915546, min_hydrophobic_side) << "TestIndexMinHydrophobicSideHelix does not give the correct results for the entry \"[PHOS]ADDAA[PHOS]D[PHOS]ADDY\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMaxHydrophobicMoment) {
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  /* peptide shorter than the window (assume w = 4), helix */
  string peptide = "A[PHOS]AA";
  double max_hmoment = RetentionFeatures::IndexMaxHydrophobicMoment(peptide, idx, 100, 4);
  EXPECT_FLOAT_EQ(0.991175806, max_hmoment) << "TestIndexMaxHydrophobicMoment does not give the correct results for the entry \"A[PHOS]AA\", window = 4, angle = 100" << endl;

  /* peptide longer than the window (assume w = 2), beta sheet */
  peptide = "A[PHOS]AD";
  max_hmoment = RetentionFeatures::IndexMaxHydrophobicMoment(peptide, idx, 180, 2);
  EXPECT_FLOAT_EQ(4.0, max_hmoment) << "TestIndexMaxHydrophobicMoment does not give the correct results for the entry \"A[PHOS]AD\", window = 2, angle = 180" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinHydrophobicMoment) {
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  /* peptide shorter than the window (assume w = 4), helix */
  string peptide = "A[PHOS]AA";
  double min_hmoment = RetentionFeatures::IndexMinHydrophobicMoment(peptide, idx, 100, 4);
  EXPECT_FLOAT_EQ(0.991175806, min_hmoment) << "TestIndexMinHydrophobicMoment does not give the correct results for the entry \"A[PHOS]AA\", window = 4, angle = 100" << endl;

  /* peptide longer than the window (assume w = 2), beta sheet */
  peptide = "A[PHOS]AD";
  min_hmoment = RetentionFeatures::IndexMinHydrophobicMoment(peptide, idx, 180, 2);
  EXPECT_FLOAT_EQ(2.0, min_hmoment) << "TestIndexMinHydrophobicMoment does not give the correct results for the entry \"A[PHOS]AD\", window = 2, angle = 180" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexSumSquaredDiff) {
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;

  string peptide = "A[PHOS]AA[PHOS]D";
  double sum_squared_diff = RetentionFeatures::IndexSumSquaredDiff(peptide, idx);
  EXPECT_FLOAT_EQ(81.0, sum_squared_diff) << "TestIndexSumSquaredDiff does not give the correct results for the entry \"A[PHOS]AA[PHOS]D\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestNumberTypeAA) {
  string temp[] = {"A", "[PHOS]D", "I"};
  set<string> type(temp, temp + 3);
  string peptide = "EEKR[PHOS]PADD[PHOS]DVVI";
  double occurences = RetentionFeatures::NumberTypeAA(peptide, type);
  EXPECT_FLOAT_EQ(3.0, occurences) << "TestNumberTypeAA does not give the correct results for the entry \"EEKR[PHOS]PADD[PHOS]DAAI\", type = A, [PHOS]D, I" << endl;
}

TEST_F(RetentionFeaturesTest, TestNumberConsecTypeAA) {
  string temp[] = {"A", "[PHOS]D", "I"};
  set<string> type(temp, temp + 3);
  string peptide = "EEKR[PHOS]PADD[PHOS]DAAVVIA";
  double occurences = RetentionFeatures::NumberConsecTypeAA(peptide, type);
  EXPECT_FLOAT_EQ(3.0, occurences) << "TestNumberConsecTypeAA does not give the correct results for the entry \"EEKR[PHOS]PADD[PHOS]DAAVVIA\", type = A, [PHOS]D, I" << endl;
}

TEST_F(RetentionFeaturesTest, TestGetExtremeRetentionAA) {
  set<string> lowest, highest;
  pair< set<string>, set<string> > set_pair;
  string exp_lowest[] = {"D", "R", "Q", "N", "K", "E"};
  string exp_highest[] = {"C", "F", "I", "L", "V"};

  set_pair = RetentionFeatures::GetExtremeRetentionAA(RetentionFeatures::k_kyte_doolittle());
  lowest = set_pair.first;
  EXPECT_EQ(6, lowest.size()) << "TestGetExtremeRetentionAA gives incorrect aa with lowest retention (length differs) " << endl;
  for (int i = 0; i < 6; ++i) {
     if (lowest.find(exp_lowest[i]) == lowest.end()) {
       ADD_FAILURE() << "TestGetExtremeRetentionAA gives incorrect aa with lowest retention " << endl;
     }
  }
  highest = set_pair.second;
  set<string>::iterator it;
  EXPECT_EQ(5, highest.size()) << "TestGetExtremeRetentionAA gives incorrect aa with highest retention (length differs) " << endl;
  for (int i = 0; i < 5; ++i) {
     if (highest.find(exp_highest[i]) == highest.end()) {
       ADD_FAILURE() << "TestGetExtremeRetentionAA gives incorrect aa with highest retention " << endl;
     }
  }
}

TEST_F(RetentionFeaturesTest, TestComputeIndexFeatures) {
  int n_features = 20;
  map<string, double> idx;
  idx["[PHOS]A"] = 5.0;
  idx["A"] = 1.0;
  idx["[PHOS]D"] = -2.0;
  idx["D"] = 3.0;
  pair< set<string>, set<string> > extreme_aa = RetentionFeatures::GetExtremeRetentionAA(idx);
  set<string> polar_aa = extreme_aa.first;
  set<string> hydrophobic_aa = extreme_aa.second;

  string peptide = "AAAAA[PHOS]A[PHOS]ADDD[PHOS]DAAAAAAAADD";
  double *features = new double[n_features];
  double expected_values[n_features];
  for(int i = 0; i < n_features; ++i)
    expected_values[i] = 0.0;
  expected_values[0] = RetentionFeatures::IndexSum(peptide, idx);
  expected_values[1] = RetentionFeatures::IndexAvg(peptide, idx);
  expected_values[2] = RetentionFeatures::IndexN(peptide, idx);
  expected_values[3] = RetentionFeatures::IndexC(peptide, idx);
  expected_values[4] = RetentionFeatures::IndexNearestNeigbour(peptide, idx, extreme_aa.first);
  expected_values[5] = RetentionFeatures::IndexMaxPartialSum(peptide, idx, 5);
  expected_values[8] = RetentionFeatures::IndexMinPartialSum(peptide, idx, 2);
  expected_values[10] = RetentionFeatures::IndexMinHydrophobicSideHelix(peptide, idx);
  expected_values[11] = RetentionFeatures::IndexMaxHydrophobicMoment(peptide, idx, 100, 11);
  expected_values[14] = RetentionFeatures::IndexMinHydrophobicMoment(peptide, idx, 180, 11);
  expected_values[15] = RetentionFeatures::IndexSumSquaredDiff(peptide, idx);
  expected_values[16] = RetentionFeatures::NumberTypeAA(peptide,  polar_aa);
  expected_values[19] = RetentionFeatures::NumberConsecTypeAA(peptide, hydrophobic_aa);
  RetentionFeatures::ComputeIndexFeatures(peptide, idx, polar_aa, hydrophobic_aa, features);
  for(int i = 0; i < n_features; ++i) {
    if (expected_values[i] != 0.0) {
     EXPECT_FLOAT_EQ(expected_values[i], *(features + i)) << "TestComputeIndexFeatures does not give the correct results for \"AAAAA[PHOS]A[PHOS]ADDD[PHOS]DAAAAAAAADD\" " << endl;
    }
  }
  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestComputeBulkinessSum) {
  string peptide = "AA[PHOS]A";

  double bulkiness_sum = RetentionFeatures::ComputeBulkinessSum(peptide, RetentionFeatures::k_bulkiness());
  EXPECT_FLOAT_EQ(34.5, bulkiness_sum) << "TestComputeBulkinessSum does not give the correct results for the entry \"AA[PHOS]A\"" << endl;
}


TEST_F(RetentionFeaturesTest, TestComputeBulkinessFeatures) {
  int n_features = 1;
  double *features = new double[n_features];
  string peptide = "AA[PHOS]A";

  RetentionFeatures::ComputeBulkinessFeatures(peptide, RetentionFeatures::k_bulkiness(), features);
  EXPECT_FLOAT_EQ(*features, RetentionFeatures::ComputeBulkinessSum(peptide, RetentionFeatures::k_bulkiness())) << "TestComputeBulkiness "
      "does not give the correct results for the entry \"AA[PHOS]A\"" << endl;

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestFillAAFeatures) {
  string temp[] = {"A", "[PHOS]D", "I", "[PHOS]A"};
  vector<string> alphabet(temp, temp + 4);
  rf.set_amino_acids_alphabet(alphabet);
  string peptide = "AAA[PHOS]D[PHOS]AII[PHOS]D";
  double expected_values[4] = {1.0, 2.0, 2.0, 3.0};
  double *aa_features = new double[4];
  aa_features = rf.FillAAFeatures(peptide, aa_features);
  for(int i = 1; i <= 4; ++i)
     EXPECT_FLOAT_EQ(expected_values[i-1], *(aa_features - i)) << "TestFillAAFeatures does not give the correct results for " << temp[3 - i] << endl;
  delete[] (aa_features - 4);
}

TEST_F(RetentionFeaturesTest, TestComputeLengthFeatures) {
  int n_features = 1;
  double *features = new double[n_features];
  string peptide = "AA[PHOS]A";

  RetentionFeatures::ComputeLengthFeatures(peptide, features);
  EXPECT_FLOAT_EQ(*features, RetentionFeatures::PeptideLength(peptide)) << "TestComputeLengthFeatures "
      "does not give the correct results for the entry \"AA[PHOS]A\"" << endl;

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestGetTotalNumberFeatures) {
  EXPECT_EQ(62, rf.GetTotalNumberFeatures());
}
TEST_F(RetentionFeaturesTest, TestComputeNoPTMFeatures) {
  rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
  int n_features = rf.GetTotalNumberFeatures();
  double *features = new double[n_features];
  string peptide = "AAAAAAAAA";

  features = rf.ComputeNoPTMFeatures(peptide, features);
  for (int i = n_features - 1; i >= 0; --i, --features) {
   if (i == n_features - 21) {
      EXPECT_FLOAT_EQ(9, *features) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
    }
   if (i == n_features - 22) {
     EXPECT_FLOAT_EQ(9, *features) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
   }
   if (i == n_features - 22) {
     EXPECT_FLOAT_EQ(9, *features) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
   }
   if (i == n_features - 26) {
     EXPECT_FLOAT_EQ(RetentionFeatures::IndexSumSquaredDiff(peptide, RetentionFeatures::k_kyte_doolittle()), *features) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
   }
  }

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestComputeRetentionFeaturesNoPtms) {
  rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
  PSMDescription psm1("AAAA", 10.0);
  PSMDescription psm2("R.YYYY.R", 11.0);
  int n_features = rf.GetTotalNumberFeatures();
  psm1.retentionFeatures = new double[n_features];
  psm2.retentionFeatures = new double[n_features];
  vector<PSMDescription> psms;
  psms.push_back(psm1);
  psms.push_back(psm2);

  rf.ComputeRetentionFeatures(psms);
  for (int i = 0; i < n_features; ++i) {
    if (i == 0) {
      EXPECT_FLOAT_EQ(RetentionFeatures::IndexSum(psm1.peptide, RetentionFeatures::k_kyte_doolittle()), psms[0].retentionFeatures[i]);
      EXPECT_FLOAT_EQ(RetentionFeatures::IndexSum(psm2.peptide.substr(2,4), RetentionFeatures::k_kyte_doolittle()), psms[1].retentionFeatures[i]);
    } if (i == 39) {
      set<string> hydrophobic_aa = RetentionFeatures::GetExtremeRetentionAA(RetentionFeatures::k_kyte_doolittle()).second;
      EXPECT_FLOAT_EQ(RetentionFeatures::NumberConsecTypeAA(psm1.peptide, hydrophobic_aa), psms[0].retentionFeatures[i]);
      EXPECT_FLOAT_EQ(RetentionFeatures::NumberConsecTypeAA(psm2.peptide.substr(2,4), hydrophobic_aa), psms[1].retentionFeatures[i]);
    } if (i == 40) {
      EXPECT_FLOAT_EQ(RetentionFeatures::ComputeBulkinessSum(psm1.peptide, RetentionFeatures::k_bulkiness()), psms[0].retentionFeatures[i]);
      EXPECT_FLOAT_EQ(RetentionFeatures::ComputeBulkinessSum(psm2.peptide.substr(2,4), RetentionFeatures::k_bulkiness()), psms[1].retentionFeatures[i]);
    } if (i == 41) {
       EXPECT_FLOAT_EQ(RetentionFeatures::PeptideLength(psm1.peptide), psms[0].retentionFeatures[i]);
       EXPECT_FLOAT_EQ(RetentionFeatures::PeptideLength(psm2.peptide.substr(2,4)), psms[1].retentionFeatures[i]);
    } if (i == 42) {
      EXPECT_FLOAT_EQ(4.0, psms[0].retentionFeatures[i]);
      EXPECT_FLOAT_EQ(0, psms[1].retentionFeatures[i]);
    } if (i == 61) {
      EXPECT_FLOAT_EQ(4.0, psms[1].retentionFeatures[i]);
      EXPECT_FLOAT_EQ(0, psms[0].retentionFeatures[i]);
    }
  }
}
