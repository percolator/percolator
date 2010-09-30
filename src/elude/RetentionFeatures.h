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
/*
 * The file stores the declaration of the class RetentionFeatures
 * This class includes information about different retention features available
 * and provides methods to calculate the retention features for a set of peptides
 */

#ifndef ELUDE_RETENTIONFEATURES_H_
#define ELUDE_RETENTIONFEATURES_H_

#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <set>
#include <utility>

/* The features are organized in groups. We will allow to switch off
 * and on only entire groups and not individual features. Here we just
 * define the index in the active_feature_groups_. As an example, if
 * active_feature_groups_[INDEX_NO_PTMS_GROUP]=1, then the group
 * NO_PTMS is switched on */
/* no ptms are present in the dataset; this is equivalent to kyte doo little, svr, peptide length, bulkiness, aa features */
#define INDEX_NO_PTMS_GROUP 0
/* phosphorylations are present in the dataset; this is equivalent to svr, peptide length, phosphorylation features, aa features  */
#define INDEX_PHOS_GROUP 1
/* only the amino acids from the alphabet are present */
#define AA_GROUP 2
/* number of feature groups */
#define NUM_FEATURE_GROUPS 3

/* forward declaration of PSMDescription */
class PSMDescription;

class RetentionFeatures {
 public:
   /* maximum number of features */
   static const int kMaxNumberFeatures;
   /* name of each feature group */
   static const std::string kGroupNames[];
   /* kyte and doolittle retention index */
   static const std::map<std::string, double> kKyteDoolittle;
   /* bulkiness as defined by Zimmerman et al */
   static const std::map<std::string, double> kBulkiness;
   /* fraction of the aa that are considered polar or hydrophobic */
   static const double kPercentageAA;

   RetentionFeatures();
   ~RetentionFeatures();

   /************ SMALL FUNCTIONS ************/
   /* get the value in the index for aa */
   static double GetIndexValue(const std::string &aa,
       const std::map<std::string, double> &index);
   /* get the amino acids in a peptide (including the modified ones) */
   static std::vector<std::string> GetAminoAcids(const std::string &peptide);
   /* get the unmodified version of an amino acid */
   static char GetUnmodifiedAA(const std::string &aa);
   /* get the total number of features */
   int GetTotalNumberFeatures() const;

   /************ INDEX FUNCTIONS ************/
   /* compute all index features for a peptide */
   static double* ComputeIndexFeatures(const std::string &peptide,
       const std::map<std::string, double> &index,
       const std::set<std::string> &polar_aa,
       const std::set<std::string> &hydrophobic_aa, double *features);
   /* get the kPercentageAA*100% AA with the lowest retention and highest retentions*/
   static std::pair< std::set<std::string>, std::set<std::string> > GetExtremeRetentionAA(
       const std::map<std::string, double> &index);
   /* calculate the number of a certain type of aa The set gives the list of such amino acids */
   static double NumberTypeAA(const std::string &peptide,
       const std::set<std::string> &amino_acid_type);
   /* calculate the number of a consecutibe aa of a certain type. The set gives the type of these such amino acids */
   static double NumberConsecTypeAA(const std::string &peptide,
       const std::set<std::string> &amino_acid_type);
   /* calculate the average hydrophobicity of an index */
   static double AvgHydrophobicityIndex(const std::map<std::string, double> &index);
   /* compute the sum of hydrophobicities of the amino acids in a peptide */
   static double IndexSum(const std::string &peptide,
       const std::map<std::string, double> &index);
   /* compute the average hydrophobicity of the aa in the peptide */
   static double IndexAvg(const std::string &peptide,
       const std::map<std::string, double> &index);
   /* calculate the hydrophobicity of the N-terminus */
   static double IndexN(const std::string &peptide, const std::map<std::string,
       double> &index);
   /* calculate the hydrophobicity of the C-terminus */
   static double IndexC(const std::string &peptide, const std::map<std::string,
       double> &index);
   /* calculate the sum of hydrophobicities of polar aa */
   static double IndexNearestNeigbour(const std::string &peptide,
       const std::map<std::string, double> &index,
       const std::set<std::string> &polar_aa);
   /* the most hydrophobic window */
   static double IndexMaxPartialSum(const std::string &peptide,
       const std::map<std::string, double> &index, const int &win);
   /* the least hydrophobic window */
   static double IndexMinPartialSum(const std::string &peptide,
       const std::map<std::string, double> &index, const int &win);
   /* calculate the most hydrophobic sides for alpha helices */
   static double IndexMaxHydrophobicSideHelix(const std::string &peptide,
       const std::map<std::string, double> &index);
   /* calculate the least hydrophobic sides for alpha helices */
   static double IndexMinHydrophobicSideHelix(const std::string &peptide,
       const std::map<std::string, double> &index);
   /* calculate the maximum value of the hydrophobic moment */
   static double IndexMaxHydrophobicMoment(const std::string &peptide,
       const std::map<std::string, double> &index, const double &angle_degrees,
       const int &win);
   /* calculate the minimum value of the hydrophobic moment */
   static double IndexMinHydrophobicMoment(const std::string &peptide,
       const std::map<std::string, double> &index, const double &angle_degrees,
       const int &win);
   /* Calculate the sum of squared differences in hydrophobicities between neighbours */
   static double IndexSumSquaredDiff(const std::string &peptide,
       const std::map<std::string, double> &index);
   /* product between hydrophobicity of n- and c- terminus */
   /* static double IndexNC(const std::string &peptide, const std::map<std::string, double> &index); */
   /* calculate the sum of hydrophobicities of neighbours of D(Aspartic Acid) and E (Glutamic acid) */
   /* static double IndexNearestNeigbourNeg(const std::string &peptide, const std::map<std::string, double> &index); */

   /************ BULKINESS FUNCTIONS ************/
   /* compute the features related to bulkiness; */
   static double* ComputeBulkinessFeatures(const std::string &peptide,
       const std::map<std::string, double> &bulkiness, double *features);
   /* compute bulkiness features; returns a pointer to the feature table */
   static double ComputeBulkinessSum(const std::string &peptide,
       const std::map<std::string, double> &bulkiness);

   /************ AMINO ACID FEATURES ************/
   /* adds a feature giving the number of each of the symbols in the alphabet found in the peptide */
   double* FillAAFeatures(const std::string &peptide, double *retention_features);

   /************ LENGTH FEATURES ************/
   /* compute the features related to length; */
   static double* ComputeLengthFeatures(const std::string &peptide, double *features);
   /* Compute the length of a peptide */
   static double PeptideLength(const std::string &peptide);

   /************* FEATURES FOR GROUPS ***************/
   /* compute the features when no ptms are present in the data */
   double* ComputeNoPTMFeatures(const std::string &peptide, double *features);
   /* compute the features when phosphorylations are present in the data */
   double* ComputePhosFeatures(const std::string &peptide, double *features);

   /************* RETENTION FEATURES FOR PSMS **************/
   /* computes the retention features for a set of peptides; return 0 if success */
   int ComputeRetentionFeatures(std::vector<PSMDescription> &psms);
   /* computes the retention features for one psm */
   int ComputeRetentionFeatures(PSMDescription &psm);

   /************* ACCESSORS AND MUTATORS **************/
   static const std::map<std::string, double>& k_kyte_doolittle()
       { return kKyteDoolittle; }
   static const std::map<std::string, double>& k_bulkiness() { return kBulkiness; }
   inline std::bitset<NUM_FEATURE_GROUPS> active_feature_groups() const
       { return active_feature_groups_; }
   inline void set_active_feature_groups(
       const std::bitset<NUM_FEATURE_GROUPS> active_groups)
       { active_feature_groups_ = active_groups; }
   inline std::map<std::string, double> svr_index() const { return svr_index_; }
   inline void set_svr_index(const std::map<std::string, double> index)
       { svr_index_ = index; }
   inline std::vector<std::string> amino_acids_alphabet() const
       { return amino_acids_alphabet_; }
   inline void set_amino_acids_alphabet(const std::vector<std::string> alphabet)
       { amino_acids_alphabet_ = alphabet; }
   static inline void set_ignore_ptms(const bool ignore_ptms) { ignore_ptms_ = ignore_ptms; }

 private:
   /* whenever a modified peptide is not identified, use the unmodified instead? */
   static bool ignore_ptms_;
   /* every bit set corresponds to an active group of features (the indices are defined at
    * the beginning of this file) */
   std::bitset<NUM_FEATURE_GROUPS> active_feature_groups_;
   /* SVR trained index */
   std::map<std::string, double> svr_index_;
   /* amino acids alphabet (can include post translationally modified peptides) */
   std::vector<std::string> amino_acids_alphabet_;
};

#endif /* ELUDE_RETENTIONFEATURES_H_ */
