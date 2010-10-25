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

#include "RetentionFeatures.h"
#include "PSMDescription.h"
#include "Globals.h"

class RetentionFeaturesTest: public ::testing::Test {
  protected:
    virtual void SetUp() {
      rf.set_ignore_ptms(true);
      Globals::getInstance()->setVerbose(1);
    }

    virtual void TearDown() {
    }

    RetentionFeatures rf;
};

TEST_F(RetentionFeaturesTest, TestGetIndexValue)
{
  /* simple amino acids */
  double val = RetentionFeatures::GetIndexValue("S",
      RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-0.8, val, 0.01)<< "GetIndexValue does not give the correct results for the entry \"S\"" << endl;

  /* modified amino acid not included in index */
  val = RetentionFeatures::GetIndexValue("Y[unimod:21]", RetentionFeatures::k_bulkiness());
  EXPECT_NEAR(18.03, val, 0.01) << "GetIndexValue does not give the correct results for the entry \"Y[unimod:21]\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["S[unimod:21]"] = 6.5;
  idx["S"] = 10.0;
  val = RetentionFeatures::GetIndexValue("S[unimod:21]", idx);
  EXPECT_NEAR(6.5, val, 0.01) << "GetIndexValue does not give the correct results for the entry \"S[unimod:21\" "
  << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestGetAminoAcids)
{
  string peptide = "A[unimod:21]NAY[unimod:21]Y[unimod:21]CC";
  string exp_peptides[7] = { "A[unimod:21]", "N", "A", "Y[unimod:21]", "Y[unimod:21]", "C", "C"};

  vector<string> res_peptides = RetentionFeatures::GetAminoAcids(peptide);
  EXPECT_EQ(7, res_peptides.size());
  for(int i = 0; i < res_peptides.size(); ++i)
    EXPECT_EQ(exp_peptides[i], res_peptides[i]);
}

TEST_F(RetentionFeaturesTest, TestIndexSum)
{
  // the index does not include the modified peptide
  string peptide = "ASA[unimod:21]Y";
  double hydrophobicity_sum = RetentionFeatures::IndexSum(peptide,
      RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(1.5, hydrophobicity_sum, 0.01)<< "IndexSum does not give the correct results for the entry \"ASA[unimod:21]AY\"" << endl;

  // modified amino acid included in index
  map<string, double> idx;
  idx["A[unimod:21]"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  hydrophobicity_sum = RetentionFeatures::IndexSum(peptide, idx);
  EXPECT_NEAR(6.7, hydrophobicity_sum, 0.01);
}

TEST_F(RetentionFeaturesTest, TestPeptideLength)
{
  // no modified aa
  double val = RetentionFeatures::PeptideLength("AASSY");
  EXPECT_NEAR(5.0, val, 0.01)<< "PeptideLength does not calculate the correct length of the peptide \"AASSY\"" << endl;

  // modified aa
  val = RetentionFeatures::PeptideLength("A[unimod:35]SY");
  EXPECT_NEAR(3.0, val, 0.01) << "PeptideLength does not calculate the correct length of the peptide \"A[unimod:35]SY\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexAvg)
{
  // the index does not include the modified peptide
  string peptide = "ASAA[unimod:21]Y";
  double hydrophobicity_avg = RetentionFeatures::IndexAvg(peptide,
      RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(3.3/5.0, hydrophobicity_avg, 0.01);

  // modified amino acid included in index
  map<string, double> idx;
  idx["C[unimod:21]"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  idx["C"] = 1.0;
  peptide = "AS[unimod:35]AC[unimod:21]Y";
  hydrophobicity_avg = RetentionFeatures::IndexAvg(peptide, idx);
  EXPECT_NEAR(7.7/5.0, hydrophobicity_avg, 0.01) << "IndexAvg does not give the correct results for the entry \"AS[unimod:35]AC[unimod:21]Y\" "
  << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbour)
{
  // get the polar aa
  set<string> polar_aa = (RetentionFeatures::GetExtremeRetentionAA(
      RetentionFeatures::k_kyte_doolittle())).first;
  // no modified aa
  string peptide = "KAYYKYCNCYYYRA";
  double sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbour(peptide,
      RetentionFeatures::k_kyte_doolittle(), polar_aa);
  EXPECT_NEAR(8.6, sum_neighbour_pos, 0.01)<< "IndexNearestNeigbour does not give the correct results for the entry \"KAYYKYCNCYYYRA\"" << endl;

  // modified aa is in the index
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["C[unimod:35]"] = -1.0;
  idx["Y[unimod:21]"] = 2.0;
  idx["A"] = 5.0;
  peptide = "KAAA[unimod:21]AC[unimod:35]A";
  polar_aa = RetentionFeatures::GetExtremeRetentionAA(idx).first;
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbour(peptide, idx, polar_aa);
  EXPECT_NEAR(10.0, sum_neighbour_pos, 0.01) << "IndexNearestNeigbour does not give the correct results for the entry \"KAAA[unimod:21]AC[unimod:35]A\"" << endl;
}

/*
 TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourPos) {
 // no modified aa
 string peptide = "KAYYKYCRCYYYRA";
 double sum_neighbour_pos =  RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(8.6, sum_neighbour_pos, 0.01) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAYYKYCRCYYYRA\"" << endl;

 // modified aa; the index does not include the modified form
 peptide = "K[unimod:21]A[unimod:21]YYKYC[unimod:21]RCYYYRA[unimod:21]";
 sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(8.6, sum_neighbour_pos, 0.01) << "IndexNearestNeigbourPos does not give the correct results for the entry \"K[unimod:21]A[unimod:21]YYKYC[unimod:21]RCYYYRA[unimod:21]\"" << endl;

 // modified aa is in the index
 map<string, double> idx;
 idx["A[unimod:21]"] = 1.0;
 idx["C[unimod:35]"] = -1.0;
 idx["Y[unimod:21]"] = 2.0;
 idx["A"] = 5.0;
 peptide = "KAAA[unimod:21]RC[unimod:35]KA";
 sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, idx);
 EXPECT_NEAR(11.0, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAAA[unimod:21]RC[unimod:35]KA\"" << endl;
 }

 TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourNeg) {
 // no modified aa
 string peptide = "DAYYDYCECYYYEA";
 double sum_neighbour_neg =  RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAYYDYCECYYYEA\"" << endl;

 // modified aa; the index does not include the modified form
 peptide = "[unimod:21]D[unimod:21]AYYDY[unimod:21]CECYYYE[unimod:21]A";
 sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"[unimod:21]D[unimod:21]AYYDY[unimod:21]CECYYYE[unimod:21]A\"" << endl;

 // modified aa is in the index
 map<string, double> idx;
 idx["[unimod:21]A"] = 1.0;
 idx["[unimod:35]C"] = -1.0;
 idx["[unimod:21]Y"] = 2.0;
 idx["A"] = 5.0;
 peptide = "DAA[unimod:21]AE[unimod:35]CDA";
 sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, idx);
 EXPECT_NEAR(11.0, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAA[unimod:21]AE[unimod:35]CDA\"" << endl;
 }*/


TEST_F(RetentionFeaturesTest, TestIndexN)
{
  // no modified aa
  string peptide = "DAYY";
  double hydrophobicity_N = RetentionFeatures::IndexN(peptide,
      RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-3.5, hydrophobicity_N, 0.01)<< "TestIndexN does not give the correct results for the entry \"DAYY\"" << endl;

  // second aa is modified but not included in the index
  peptide = "YD[unimod:21]AYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-1.3, hydrophobicity_N, 0.01) << "TestIndexN does not give the correct results for the entry \"YD[unimod:21]AYY\"" << endl;

  // N term is modified but not included in the index
  peptide = "D[unimod:21]AYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-3.5, hydrophobicity_N, 0.01) << "TestIndexN does not give the correct results for the entry \"D[unimod:21]AYY\"" << endl;

  // N term is modified and included in the index
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  peptide = "A[unimod:21]E";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, idx);
  EXPECT_NEAR(1.0, hydrophobicity_N, 0.01) << "TestIndexN does not give the correct results for the entry \"A[unimod:21]E\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexC)
{
  // no modified aa
  string peptide = "DAYY";
  double hydrophobicity_C = RetentionFeatures::IndexC(peptide,
      RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-1.3, hydrophobicity_C, 0.01)<< "TestIndexC does not give the correct results for the entry \"DAYY\"" << endl;

  // one but last is modified but not included in the index
  peptide = "YD[unimod:21]AYYD[unimod:21]D";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-3.5, hydrophobicity_C, 0.01) << "TestIndexC does not give the correct results for the entry \"YD[unimod:21]AYYD[unimod:21]D\"" << endl;

  // C term is modified but not included in the index
  peptide = "D[unimod:21]AYD[unimod:21]";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_NEAR(-3.5, hydrophobicity_C, 0.01) << "TestIndexC does not give the correct results for the entry \"D[unimod:21]AYD[unimod:21]\"" << endl;

  // C term is modified and included in the index
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  peptide = "EA[unimod:21]";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, idx);
  EXPECT_NEAR(1.0, hydrophobicity_C, 0.01) << "TestIndexC does not give the correct results for the entry \"EA[unimod:21]\"" << endl;
}

/*
 TEST_F(RetentionFeaturesTest, TestIndexNC) {
 // no modified aa
 string peptide = "DAYY";
 double hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(0.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"DAYY\"" << endl;

 // second and last aa are modified but not included in the index
 peptide = "I[unimod:21]DAYY[unimod:21]F";
 hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"I[unimod:21]DAYY[unimod:21]F\"" << endl;

 // both c and n terminus are modified but not included in the index
 peptide = "[unimod:21]IAY[unimod:21]F";
 hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
 EXPECT_NEAR(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[unimod:21]IAY[unimod:21]F\"" << endl;

 // C term and N term are modified and included in the index
 map<string, double> idx;
 idx["[unimod:21]A"] = 1.0;
 idx["A"] = 5.0;
 peptide = "[unimod:21]A";
 hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, idx);
 EXPECT_NEAR(1.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[unimod:21]A\"" << endl;
 }*/

TEST_F(RetentionFeaturesTest, TestIndexMaxPartialSum)
{
  // no modified aa
  string peptide = "DDDFFFDDDIIDDD";
  double max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide,
      RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_NEAR(9.0, max_partial_sum, 0.01)<< "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  // modified aa but not in the index
  peptide = "DDDAA[unimod:21]ADDDII[unimod:21]WDD";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_NEAR(8.1, max_partial_sum, 0.01) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDAA[unimod:21]ADDDII[unimod:21]WDD\"" << endl;

  // modified aa that are in the index
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = 2.0;
  idx["D"] = 3.0;
  peptide = "A[unimod:21]ADDD[unimod:21]";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, idx, 3);
  EXPECT_NEAR(11.0, max_partial_sum, 0.01) << "TestIndexMaxPartialSum does not give the correct results for the entry \"A[unimod:21]ADDD[unimod:21]\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinPartialSum)
{
  // no modified aa
  string peptide = "DDDFFFDDDIIDDD";
  double min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide,
      RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_NEAR(-7.0, min_partial_sum, 0.01)<< "TestIndexMinPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  // modified aa but not in the index
  peptide = "D[unimod:21]D[unimod:21]DADD";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_NEAR(-10.5, min_partial_sum, 0.01) << "TestIndexMinPartialSum does not give the correct results for the entry \"D[unimod:21]D[unimod:21]DADD\"" << endl;

  // modified aa that are in the index
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = 2.0;
  idx["D"] = 3.0;
  peptide = "A[unimod:21]ADDD[unimod:21]";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, idx, 3);
  EXPECT_NEAR(8.0, min_partial_sum, 0.01) << "TestIndexMinPartialSum does not give the correct results for the entry \"A[unimod:21]ADDD[unimod:21]\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestAvgHydrophobicityIndex)
{
  // no modified aa
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = 2.0;
  idx["D"] = 3.0;
  double avg_hydrophobicity = RetentionFeatures::AvgHydrophobicityIndex(idx);
  EXPECT_NEAR(2.75, avg_hydrophobicity, 0.01) << "TestAvgHydrophobicityIndex does not give the correct result" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMaxHydrophobicSideHelix)
{
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  // peptide smaller than 9 aa
  string peptide = "DDDDDD";
  double max_hydrophobic_side =
      RetentionFeatures::IndexMaxHydrophobicSideHelix(peptide, idx);
  EXPECT_NEAR(7.841237327, max_hydrophobic_side, 0.01) << "TestIndexMaxHydrophobicSideHelix does not give the correct results for the entry \"DDDDDD\"" << endl;

  // normal peptides (with modified aa)
  peptide = "A[unimod:21]DDAAD[unimod:21]A[unimod:21]DDY";
  max_hydrophobic_side = RetentionFeatures::IndexMaxHydrophobicSideHelix(peptide, idx);
  EXPECT_NEAR(11.064177772, max_hydrophobic_side, 0.01) << "TestIndexMaxHydrophobicSideHelix does not give the correct results for the entry \"A[unimod:21]DDAAD[unimod:21]A[unimod:21]DDY\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinHydrophobicSideHelix)
{
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  // peptide smaller than 9 aa
  string peptide = "DDDDDD";
  double min_hydrophobic_side =
      RetentionFeatures::IndexMinHydrophobicSideHelix(peptide, idx);
  EXPECT_NEAR(7.841237327, min_hydrophobic_side, 0.01)<< "TestIndexMinHydrophobicSideHelix does not give the correct results for the entry \"DDDDDD\"" << endl;

  // normal peptides (with modified aa)
  peptide = "A[unimod:21]DDAAD[unimod:21]A[unimod:21]DDY";
  min_hydrophobic_side = RetentionFeatures::IndexMinHydrophobicSideHelix(peptide, idx);
  EXPECT_NEAR(6.438915546, min_hydrophobic_side, 0.01) << "TestIndexMinHydrophobicSideHelix does not give the correct results for the entry \"A[unimod:21]DDAAD[unimod:21]A[unimod:21]DDY\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMaxHydrophobicMoment)
{
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  // peptide shorter than the window (assume w = 4), helix
  string peptide = "AA[unimod:21]A";
  double max_hmoment = RetentionFeatures::IndexMaxHydrophobicMoment(peptide,
      idx, 100, 4);
  EXPECT_NEAR(0.991175806, max_hmoment, 0.01)<< "TestIndexMaxHydrophobicMoment does not give the correct results for the entry \"AA[unimod:21]A\", window = 4, angle = 100" << endl;

  // peptide longer than the window (assume w = 2), beta sheet
  peptide = "AA[unimod:21]D";
  max_hmoment = RetentionFeatures::IndexMaxHydrophobicMoment(peptide, idx, 180, 2);
  EXPECT_NEAR(4.0, max_hmoment, 0.01) << "TestIndexMaxHydrophobicMoment does not give the correct results for the entry \"AA[unimod:21]D\", window = 2, angle = 180" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinHydrophobicMoment)
{
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;
  idx["Y"] = 4.1;

  // peptide shorter than the window (assume w = 4), helix
  string peptide = "AAA[unimod:21]";
  double min_hmoment = RetentionFeatures::IndexMinHydrophobicMoment(peptide,
      idx, 100, 4);
  EXPECT_NEAR(0.991175806, min_hmoment, 0.01)<< "TestIndexMinHydrophobicMoment does not give the correct results for the entry \"A[unimod:21]AA\", window = 4, angle = 100" << endl;

  // peptide longer than the window (assume w = 2), beta sheet
  peptide = "A[unimod:21]AD";
  min_hmoment = RetentionFeatures::IndexMinHydrophobicMoment(peptide, idx, 180, 2);
  EXPECT_NEAR(2.0, min_hmoment, 0.01) << "TestIndexMinHydrophobicMoment does not give the correct results for the entry \"A[unimod:21]AD\", window = 2, angle = 180" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexSumSquaredDiff)
{
  map<string, double> idx;
  idx["A[unimod:21]"] = 1.0;
  idx["A"] = 5.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;

  string peptide = "AA[unimod:21]AD[unimod:21]";
  double sum_squared_diff =
      RetentionFeatures::IndexSumSquaredDiff(peptide, idx);
  EXPECT_NEAR(81.0, sum_squared_diff, 0.01)<< "TestIndexSumSquaredDiff does not give the correct results for the entry \"AA[unimod:21]AD[unimod:21]\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestNumberTypeAA)
{
  string temp[] = { "A", "D[unimod:21]", "I" };
  set<string> type(temp, temp + 3);
  string peptide = "EEKRP[unimod:21]ADDD[unimod:21]VVI";
  double occurences = RetentionFeatures::NumberTypeAA(peptide, type);
  EXPECT_NEAR(3.0, occurences, 0.0001)<< "TestNumberTypeAA does not give the correct results for the entry \"EEKRP[unimod:21]ADDD[unimod:21]VVI\", type = A, D[unimod:21], I" << endl;
}

TEST_F(RetentionFeaturesTest, TestNumberConsecTypeAA)
{
  string temp[] = { "A", "D[unimod:21]", "I" };
  set<string> type(temp, temp + 3);
  string peptide = "EEKRP[unimod:21]ADDD[unimod:21]AAVVIA";
  double occurences = RetentionFeatures::NumberConsecTypeAA(peptide, type);
  EXPECT_NEAR(3.0, occurences, 0.001)<< "TestNumberConsecTypeAA does not give the correct results for the entry \"EEKRP[unimod:21]ADDD[unimod:21]AAVVIA\", type = A, D[unimod:21], I" << endl;
}

TEST_F(RetentionFeaturesTest, TestGetExtremeRetentionAA)
{
  set<string> lowest, highest;
  pair<set<string> , set<string> > set_pair;
  string exp_lowest[] = { "D", "R", "Q", "N", "K", "E" };
  string exp_highest[] = { "C", "F", "I", "L", "V" };

  set_pair = RetentionFeatures::GetExtremeRetentionAA(
      RetentionFeatures::k_kyte_doolittle());
  lowest = set_pair.first;
  EXPECT_EQ(6, lowest.size())<< "TestGetExtremeRetentionAA gives incorrect aa with lowest retention (length differs) " << endl;
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

TEST_F(RetentionFeaturesTest, TestComputeIndexFeatures)
{
  int n_features = 20;
  map<string, double> idx;
  idx["A[unimod:21]"] = 5.0;
  idx["A"] = 1.0;
  idx["D[unimod:21]"] = -2.0;
  idx["D"] = 3.0;
  pair<set<string> , set<string> > extreme_aa =
      RetentionFeatures::GetExtremeRetentionAA(idx);
  set<string> polar_aa = extreme_aa.first;
  set<string> hydrophobic_aa = extreme_aa.second;

  string peptide = "AAAAAA[unimod:21]A[unimod:21]DDDD[unimod:21]AAAAAAAADD";
  double *features = new double[n_features];
  double expected_values[n_features];
  for (int i = 0; i < n_features; ++i)
    expected_values[i] = 0.0;
  expected_values[0] = RetentionFeatures::IndexSum(peptide, idx);
  expected_values[1] = RetentionFeatures::IndexAvg(peptide, idx);
  expected_values[2] = RetentionFeatures::IndexN(peptide, idx);
  expected_values[3] = RetentionFeatures::IndexC(peptide, idx);
  expected_values[4] = RetentionFeatures::IndexNearestNeigbour(peptide, idx,
      extreme_aa.first);
  expected_values[5] = RetentionFeatures::IndexMaxPartialSum(peptide, idx, 5);
  expected_values[8] = RetentionFeatures::IndexMinPartialSum(peptide, idx, 2);
  expected_values[10] = RetentionFeatures::IndexMinHydrophobicSideHelix(
      peptide, idx);
  expected_values[11] = RetentionFeatures::IndexMaxHydrophobicMoment(peptide,
      idx, 100, 11);
  expected_values[14] = RetentionFeatures::IndexMinHydrophobicMoment(peptide,
      idx, 180, 11);
  expected_values[15] = RetentionFeatures::IndexSumSquaredDiff(peptide, idx);
  expected_values[16] = RetentionFeatures::NumberTypeAA(peptide, polar_aa);
  expected_values[19] = RetentionFeatures::NumberConsecTypeAA(peptide,
      hydrophobic_aa);
  RetentionFeatures::ComputeIndexFeatures(peptide, idx, polar_aa,
      hydrophobic_aa, features);
  for (int i = 0; i < n_features; ++i) {
    if (expected_values[i] != 0.0) {
      EXPECT_NEAR(expected_values[i], *(features + i), 0.01)<< "TestComputeIndexFeatures does not give the correct results for \"AAAAAA[unimod:21]A[unimod:21]DDDD[unimod:21]AAAAAAAADD\" " << endl;
    }
  }
  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestComputeBulkinessSum)
{
  string peptide = "AAA[unimod:21]";

  double bulkiness_sum = RetentionFeatures::ComputeBulkinessSum(peptide,
      RetentionFeatures::k_bulkiness());
  EXPECT_NEAR(34.5, bulkiness_sum, 0.01)<< "TestComputeBulkinessSum does not give the correct results for the entry \"AAA[unimod:21]\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestComputeBulkinessFeatures)
{
  int n_features = 1;
  double *features = new double[n_features];
  string peptide = "AA[unimod:21]A";

  RetentionFeatures::ComputeBulkinessFeatures(peptide,
      RetentionFeatures::k_bulkiness(), features);
  EXPECT_NEAR(*features, RetentionFeatures::ComputeBulkinessSum(peptide, RetentionFeatures::k_bulkiness()), 0.01)<< "TestComputeBulkiness "
  "does not give the correct results for the entry \"AA[unimod:21]A\"" << endl;

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestFillAAFeatures)
{
  string temp[] = { "A", "D[unimod:21]", "I", "A[unimod:21]" };
  vector<string> alphabet(temp, temp + 4);
  rf.set_amino_acids_alphabet(alphabet);
  string peptide = "AAAD[unimod:21]A[unimod:21]IID[unimod:21]";
  double expected_values[4] = { 1.0, 2.0, 2.0, 3.0 };
  double *aa_features = new double[4];
  aa_features = rf.FillAAFeatures(peptide, aa_features);
  for (int i = 1; i <= 4; ++i)
    EXPECT_NEAR(expected_values[i-1], *(aa_features - i), 0.01)<< "TestFillAAFeatures does not give the correct results for " << temp[3 - i] << endl;
    delete[] (aa_features - 4);
}

TEST_F(RetentionFeaturesTest, TestComputeLengthFeatures)
{
  int n_features = 1;
  double *features = new double[n_features];
  string peptide = "AA[unimod:21]A";

  RetentionFeatures::ComputeLengthFeatures(peptide, features);
  EXPECT_NEAR(*features, RetentionFeatures::PeptideLength(peptide), 0.01)<< "TestComputeLengthFeatures "
      "does not give the correct results for the entry \"AA[unimod:21]A\"" << endl;

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestGetTotalNumberFeatures)
{
  EXPECT_EQ(62, rf.GetTotalNumberFeatures());
}

TEST_F(RetentionFeaturesTest, TestComputeNoPTMFeatures)
{
  rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
  int n_features = rf.GetTotalNumberFeatures();
  double *features = new double[n_features];
  string peptide = "AAAAAAAAA";

  features = rf.ComputeNoPTMFeatures(peptide, features);
  for (int i = n_features - 1; i >= 0; --i, --features) {
    if (i == n_features - 21) {
      EXPECT_NEAR(9, *features, 0.01)<< "TestComputeNoPTMFeatures does not give the correct result" << endl;
    }
    if (i == n_features - 22) {
      EXPECT_NEAR(9, *features, 0.01) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
    }
    if (i == n_features - 22) {
      EXPECT_NEAR(9, *features, 0.01) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
    }
    if (i == n_features - 26) {
      EXPECT_NEAR(RetentionFeatures::IndexSumSquaredDiff(peptide, RetentionFeatures::k_kyte_doolittle()), *features, 0.01) << "TestComputeNoPTMFeatures does not give the correct result" << endl;
    }
  }

  delete[] features;
}

TEST_F(RetentionFeaturesTest, TestComputeRetentionFeaturesNoPtms)
{
  rf.set_svr_index(RetentionFeatures::k_kyte_doolittle());
  PSMDescription psm1("AAAA[unimod:21]", 10.0);
  PSMDescription psm2("R.YY[unimod:21]YY.R", 11.0);
  int n_features = rf.GetTotalNumberFeatures();
  psm1.retentionFeatures = new double[n_features];
  psm2.retentionFeatures = new double[n_features];
  vector<PSMDescription> psms;
  psms.push_back(psm1);
  psms.push_back(psm2);

  rf.ComputeRetentionFeatures(psms);
  for (int i = 0; i < n_features; ++i) {
    if (i == 0) {
      EXPECT_NEAR(RetentionFeatures::IndexSum(psm1.peptide, RetentionFeatures::k_kyte_doolittle()), psms[0].retentionFeatures[i], 0.01) << " i = 0";
      EXPECT_NEAR(RetentionFeatures::IndexSum(psm2.peptide.substr(2,15), RetentionFeatures::k_kyte_doolittle()), psms[1].retentionFeatures[i], 0.01)  << " i = 0";
    } if (i == 39) {
      set<string> hydrophobic_aa = RetentionFeatures::GetExtremeRetentionAA(RetentionFeatures::k_kyte_doolittle()).second;
      EXPECT_NEAR(RetentionFeatures::NumberConsecTypeAA(psm1.peptide, hydrophobic_aa), psms[0].retentionFeatures[i], 0.01)  << " i = 39";
      EXPECT_NEAR(RetentionFeatures::NumberConsecTypeAA(psm2.peptide.substr(2,15), hydrophobic_aa), psms[1].retentionFeatures[i], 0.01) << " i = 39";
    }if (i == 40) {
      EXPECT_NEAR(RetentionFeatures::ComputeBulkinessSum(psm1.peptide, RetentionFeatures::k_bulkiness()), psms[0].retentionFeatures[i], 0.01) << " i = 40";
      EXPECT_NEAR(RetentionFeatures::ComputeBulkinessSum(psm2.peptide.substr(2,15), RetentionFeatures::k_bulkiness()), psms[1].retentionFeatures[i], 0.01) << " i = 40";
    }if (i == 41) {
      EXPECT_NEAR(RetentionFeatures::PeptideLength(psm1.peptide), psms[0].retentionFeatures[i], 0.01) << " i = 41";
      EXPECT_NEAR(RetentionFeatures::PeptideLength(psm2.peptide.substr(2,15)), psms[1].retentionFeatures[i], 0.01) << " i = 41";
    }if (i == 42) {
      EXPECT_NEAR(4.0, psms[0].retentionFeatures[i], 0.01) << " i = 42";
      EXPECT_FLOAT_EQ(0, psms[1].retentionFeatures[i]) << " i = 42";
    }if (i == 61) {
      EXPECT_NEAR(4.0, psms[1].retentionFeatures[i], 0.01) << " i = 61";
      EXPECT_FLOAT_EQ(0, psms[0].retentionFeatures[i]) << " i = 61";
    }
  }
}
