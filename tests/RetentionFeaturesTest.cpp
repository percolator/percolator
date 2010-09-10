
#include <gtest/gtest.h>
#include <vector>
#include <map>
#include <set>

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
  EXPECT_DOUBLE_EQ(3.3, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  hydrophobicity_sum =  RetentionFeatures::IndexSum(peptide, idx);
  EXPECT_DOUBLE_EQ(7.7, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\" "
                                            << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestPeptideLength) {
  /* no modified aa */
  double val = RetentionFeatures::PeptideLength("AASSY");
  EXPECT_DOUBLE_EQ(5.0, val) << "PeptideLength does not calculate the correct length of the peptide \"AASSY\"" << endl;

  /* modified aa */
  val = RetentionFeatures::PeptideLength("[PHOS]AS[GLYC]Y");
  EXPECT_DOUBLE_EQ(3.0, val) << "PeptideLength does not calculate the correct length of the peptide \"[PHOS]A[GLYC]Y\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexAvg) {
  /* the index does not include the modified peptide */
  string peptide = "ASA[PHOS]AY";
  double hydrophobicity_avg =  RetentionFeatures::IndexAvg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(3.3/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  peptide = "A[GLYC]SA[PHOS]AY";
  hydrophobicity_avg =  RetentionFeatures::IndexAvg(peptide, idx);
  EXPECT_DOUBLE_EQ(7.7/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"A[GLYC]SA[PHOS]AY\" "
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

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourPos) {
  /* no modified aa */
  string peptide = "KAYYKYCRCYYYRA";
  double sum_neighbour_pos =  RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAYYKYCRCYYYRA\"" << endl;

  /* modified aa; the index does not include the modified form */
  peptide = "[PHOS]K[PHOS]AYYKY[PHOS]CRCYYYR[PHOS]A";
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"[PHOS]K[PHOS]AYYKY[PHOS]CRCYYYR[PHOS]A\"" << endl;

  /* modified aa is in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["[GLYC]C"] = -1.0;
  idx["[PHOS]Y"] = 2.0;
  idx["A"] = 5.0;
  peptide = "KAA[PHOS]AR[GLYC]CKA";
  sum_neighbour_pos = RetentionFeatures::IndexNearestNeigbourPos(peptide, idx);
  EXPECT_DOUBLE_EQ(11.0, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAA[PHOS]AR[GLYC]CKA\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourNeg) {
  /* no modified aa */
  string peptide = "DAYYDYCECYYYEA";
  double sum_neighbour_neg =  RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAYYDYCECYYYEA\"" << endl;

  /* modified aa; the index does not include the modified form */
  peptide = "[PHOS]D[PHOS]AYYDY[PHOS]CECYYYE[PHOS]A";
  sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"[PHOS]D[PHOS]AYYDY[PHOS]CECYYYE[PHOS]A\"" << endl;

  /* modified aa is in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["[GLYC]C"] = -1.0;
  idx["[PHOS]Y"] = 2.0;
  idx["A"] = 5.0;
  peptide = "DAA[PHOS]AE[GLYC]CDA";
  sum_neighbour_neg = RetentionFeatures::IndexNearestNeigbourNeg(peptide, idx);
  EXPECT_DOUBLE_EQ(11.0, sum_neighbour_neg) << "IndexNearestNeigbourNeg does not give the correct results for the entry \"DAA[PHOS]AE[GLYC]CDA\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexN) {
  /* no modified aa */
  string peptide = "DAYY";
  double hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-3.5, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"DAYY\"" << endl;

  /* second aa is modified but not included in the index */
  peptide = "Y[PHOS]DAYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-1.3, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"Y[PHOS]DAYY\"" << endl;

  /* N term is modified but not included in the index */
  peptide = "[PHOS]DAYY";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-3.5, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"[PHOS]DAYY\"" << endl;

  /* N term is modified and included in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "[PHOS]AE";
  hydrophobicity_N = RetentionFeatures::IndexN(peptide, idx);
  EXPECT_DOUBLE_EQ(1.0, hydrophobicity_N) << "TestIndexN does not give the correct results for the entry \"[PHOS]AE\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexC) {
  /* no modified aa */
  string peptide = "DAYY";
  double hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-1.3, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"DAYY\"" << endl;

  /* one but last is modified but not included in the index */
  peptide = "Y[PHOS]DAYY[PHOS]DD";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-3.5, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"Y[PHOS]DAYY[PHOS]DD\"" << endl;

  /* C term is modified but not included in the index */
  peptide = "[PHOS]DAY[PHOS]D";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(-3.5, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"[PHOS]DAY[PHOS]D\"" << endl;

  /* C term is modified and included in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "E[PHOS]A";
  hydrophobicity_C = RetentionFeatures::IndexC(peptide, idx);
  EXPECT_DOUBLE_EQ(1.0, hydrophobicity_C) << "TestIndexC does not give the correct results for the entry \"E[PHOS]A\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexNC) {
  /* no modified aa */
  string peptide = "DAYY";
  double hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(0.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"DAYY\"" << endl;

  /* second and last aa are modified but not included in the index */
  peptide = "I[PHOS]DAYY[PHOS]F";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"I[PHOS]DAYY[PHOS]F\"" << endl;

  /* both c and n terminus are modified but not included in the index */
  peptide = "[PHOS]IAY[PHOS]F";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(12.6, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[PHOS]IAY[PHOS]F\"" << endl;

  /* C term and N term are modified and included in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  peptide = "[PHOS]A";
  hydrophobicity_NC = RetentionFeatures::IndexNC(peptide, idx);
  EXPECT_DOUBLE_EQ(1.0, hydrophobicity_NC) << "TestIndexNC does not give the correct results for the entry \"[PHOS]A\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMaxPartialSum) {
  /* no modified aa */
  string peptide = "DDDFFFDDDIIDDD";
  double max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_DOUBLE_EQ(9.0, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  /* modified aa but not in the index */
  peptide = "DDDA[PHOS]AADDDI[PHOS]IWDD";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_DOUBLE_EQ(8.1, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;

  /* modified aa that are in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  peptide = "[PHOS]AADD[PHOS]D";
  max_partial_sum = RetentionFeatures::IndexMaxPartialSum(peptide, idx, 3);
  EXPECT_DOUBLE_EQ(11.0, max_partial_sum) << "TestIndexMaxPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexMinPartialSum) {
  /* no modified aa */
  string peptide = "DDDFFFDDDIIDDD";
  double min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 2);
  EXPECT_DOUBLE_EQ(-7.0, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDFFFDDDIIDDD\"" << endl;

  /* modified aa but not in the index */
  peptide = "[PHOS]D[PHOS]DDADD";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, RetentionFeatures::k_kyte_doolittle(), 3);
  EXPECT_DOUBLE_EQ(-10.5, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;

  /* modified aa that are in the index */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  peptide = "[PHOS]AADD[PHOS]D";
  min_partial_sum = RetentionFeatures::IndexMinPartialSum(peptide, idx, 3);
  EXPECT_DOUBLE_EQ(8.0, min_partial_sum) << "TestIndexMinPartialSum does not give the correct results for the entry \"DDDA[PHOS]AADDDI[PHOS]IWDD\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestAvgHydrophobicityIndex) {
  /* no modified aa */
  map<string, double> idx;
  idx["[PHOS]A"] = 1.0;
  idx["A"] = 5.0;
  idx["[PHOS]D"] = 2.0;
  idx["D"] = 3.0;
  double avg_hydrophobicity = RetentionFeatures:: AvgHydrophobicityIndex(idx);
  EXPECT_DOUBLE_EQ(2.75, avg_hydrophobicity) << "TestAvgHydrophobicityIndex does not give the correct result" << endl;
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


FillAAFeatures(const string &peptide, double *retention_features)
