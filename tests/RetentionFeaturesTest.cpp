
#include <gtest/gtest.h>
#include <vector>
#include <map>

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
  double val = rf.GetIndexValue("S", RetentionFeatures::k_kyte_doolittle());
  EXPECT_EQ(-0.8, val) << "GetIndexValue does not give the correct results for the entry \"S\"" << endl;

  /* modified amino acid not included in index */
  val = rf.GetIndexValue("[PHOS]Y", RetentionFeatures::k_bulkiness());
  EXPECT_EQ(18.03, val) << "GetIndexValue does not give the correct results for the entry \"[PHOS]Y\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 6.5;
  idx["A"] = 10.0;
  val = rf.GetIndexValue("[PHOS]A", idx);
  EXPECT_EQ(6.5, val) << "GetIndexValue does not give the correct results for the entry \"[PHOS]A\" "
                      << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexSum) {
  /* the index does not include the modified peptide */
  string peptide = "ASA[PHOS]AY";
  double hydrophobicity_sum = rf.IndexSum(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(3.3, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  hydrophobicity_sum = rf.IndexSum(peptide, idx);
  EXPECT_DOUBLE_EQ(7.7, hydrophobicity_sum) << "IndexSum does not give the correct results for the entry \"AA[PHOSPH]AY\" "
                                            << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestPeptideLength) {
  /* no modified aa */
  double val = rf.PeptideLength("AASSY");
  EXPECT_DOUBLE_EQ(5.0, val) << "PeptideLength does not calculate the correct length of the peptide \"AASSY\"" << endl;

  /* modified aa */
  val = rf.PeptideLength("[PHOS]AS[GLYC]Y");
  EXPECT_DOUBLE_EQ(3.0, val) << "PeptideLength does not calculate the correct length of the peptide \"[PHOS]A[GLYC]Y\"" << endl;
}

TEST_F(RetentionFeaturesTest, TestIndexAvg) {
  /* the index does not include the modified peptide */
  string peptide = "ASA[PHOS]AY";
  double hydrophobicity_avg = rf.IndexAvg(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(3.3/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"AA[PHOSPH]AY\"" << endl;

  /* modified amino acid included in index */
  map<string, double> idx;
  idx["[PHOS]A"] = 0.7;
  idx["A"] = 1.0;
  idx["Y"] = 2.0;
  idx["S"] = 3.0;
  peptide = "A[GLYC]SA[PHOS]AY";
  hydrophobicity_avg = rf.IndexAvg(peptide, idx);
  EXPECT_DOUBLE_EQ(7.7/5.0, hydrophobicity_avg) << "IndexAvg does not give the correct results for the entry \"A[GLYC]SA[PHOS]AY\" "
                                            << "(index includes modified peptide)" << endl;
}

TEST_F(RetentionFeaturesTest, TestGetAminoAcids) {
  string peptide = "[PHOS]AAA[PHOS]Y[PHOS]Y";
  string exp_peptides[5] = {"[PHOS]A", "A", "A", "[PHOS]Y", "[PHOS]Y"};

  vector<string> res_peptides = rf.GetAminoAcids(peptide);
  EXPECT_EQ(res_peptides.size(), 5);
  for(int i = 0; i < res_peptides.size(); ++i)
    EXPECT_EQ(exp_peptides[i], res_peptides[i]);
}

TEST_F(RetentionFeaturesTest, TestIndexNearestNeigbourPos) {
  /* no modified aa */
  string peptide = "KAYYKYCRCYYYRA";
  double sum_neighbour_pos = rf.IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAYYKYCRCYYYRA\"" << endl;

  /* modified aa; the index does not include the modified form */
  //peptide = "[PHOS]K[PHOS]AYYKY[PHOS]CRCYYYR[PHOS]A";
  sum_neighbour_pos = rf.IndexNearestNeigbourPos(peptide, RetentionFeatures::k_kyte_doolittle());
  EXPECT_DOUBLE_EQ(8.6, sum_neighbour_pos) << "IndexNearestNeigbourPos does not give the correct results for the entry \"KAYYKYCRCYYYRA\"" << endl;

}
