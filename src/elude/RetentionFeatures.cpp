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
 * The file includes definitions of variables and methods in the class RetentionFeatures
 */
#include <math.h>
#include <stdio.h>

#include <iostream>
#include <algorithm>

#include "boost/assign.hpp"
#include "RetentionFeatures.h"
#include "Globals.h"
#include "PSMDescription.h"

using namespace boost::assign;

/* maximum number of features */
const int kMaxNumberFeatures = 100;

/* define the Kyte Doolittle index */
const map<string, double> RetentionFeatures::kKyteDoolittle = map_list_of ("A", 1.8) ("C", 2.5) ("D", -3.5) ("E", -3.5) ("F", 2.8)
                                                                          ("G", -0.4) ("H", -3.2) ("I", 4.5) ("K", -3.9) ("L", 3.8)
                                                                          ("M", 1.9) ("N", -3.5) ("P", -1.6) ("Q", -3.5) ("R", -4.5)
                                                                          ("S", -0.8) ("T", -0.7) ("V", 4.2) ("W", -0.9) ("Y", -1.3);
/* define the bulkiness */
const map<string, double> RetentionFeatures::kBulkiness = map_list_of ("A", 11.5) ("C", 13.46) ("D", 11.68) ("E", 13.57) ("F", 19.80)
                                                                      ("G", 3.40) ("H", 13.69) ("I", 21.40) ("K", 15.71) ("L", 21.40)
                                                                      ("M", 16.25) ("N", 12.82) ("P", 17.43) ("Q", 14.45) ("R", 14.28)
                                                                      ("S", 9.47) ("T", 15.77) ("V", 21.57) ("W", 21.67) ("Y", 18.03);

/* percentage of the amino acids from an index that are polar or hydrophobic */
const double RetentionFeatures::kPercentageAA = 0.25;

/* name of each feature group */
const string RetentionFeatures::kGroupNames[] = {"No posttranslationally modified peptides", "Phosphorylations"};

RetentionFeatures::RetentionFeatures() {
  string aa_alphabet[] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};
  amino_acids_alphabet_.assign(aa_alphabet, aa_alphabet + 20);
  // by default we have no ptms
  active_feature_groups_.set(INDEX_NO_PTMS_GROUP);
}

RetentionFeatures::~RetentionFeatures() {
}

/**************************** SMALL FUNCTIONS **************************************/
/* Return the index value of an modified amino acid; if it is not included in the index,
 * then the value of the unmodified one is returned and an warning is returned */
double RetentionFeatures::GetIndexValue(const string &aa, const map<string, double> &index) {
  map<string, double>::const_iterator index_value = index.find(aa);
  string unmodified_aa;

  if (index_value != index.end()) {
    return index_value->second;
  } else {
    unmodified_aa = aa[aa.size() - 1];
    index_value = index.find(unmodified_aa);
    if (index_value != index.end()) {
      if (VERB >= 4) {
        cerr << "Warning: Could not find the index value for " << aa << ". Use index value for " << unmodified_aa << " instead" << endl;
      }
      return index_value->second;
    } else {
      if (VERB >= 1) {
        cerr << "Error: Could not find the index value for " << aa << endl;
        cerr << "Execution aborted" << endl;
      }
      exit(1);
    }
  }
}

/* get the amino acids in a peptide */
vector<string> RetentionFeatures::GetAminoAcids(const string &peptide) {
  vector<string> amino_acids;
  int end_position, i = 0, len = peptide.size();
  string aa;

  while (i < len) {
     aa = peptide.at(i);
     if (aa == "[") {
       end_position = peptide.find("]", i);
       aa = peptide.substr(i, end_position - i + 2);
     } else if (aa == "-") {
       i = i + 1;
       continue;
     }
     amino_acids.push_back(aa);
     i += aa.size();
   }
  return amino_acids;
}

/* get the unmodified amino acid */
char RetentionFeatures::GetUnmodifiedAA(const string &aa) {
  return aa[aa.size() - 1];
}

/* get the total number of features */
int RetentionFeatures::GetTotalNumberFeatures() const {
  int number_index_features, number_length_features, number_aa_features;
  if (active_feature_groups_.test(INDEX_NO_PTMS_GROUP)) {
    int number_bulkiness_features = 1;
    number_index_features = 20;
    number_length_features = 1;
    number_aa_features = 20;
    return 2 * number_index_features +  number_length_features + number_bulkiness_features +  number_aa_features;
  } else {
    //[TO DO: modify when more features are added]
    number_index_features = 20;
    number_length_features = 1;
    number_aa_features = amino_acids_alphabet_.size();

    return number_index_features +  number_length_features +  number_aa_features;
  }

}

/**************************** INDEX FUNCTIONS **************************************/
/* fill all the features of an index; return a pointer to the next element in the feature table */
double* RetentionFeatures::ComputeIndexFeatures(const string &peptide, const map<string, double> &index, const set<string> &polar_aa,
                           const set<string> &hydrophobic_aa, double *features) {
  *(features++) = IndexSum(peptide, index);
  *(features++) = IndexAvg(peptide, index);
  *(features++) = IndexN(peptide, index);
  *(features++) = IndexC(peptide, index);
  *(features++) = IndexNearestNeigbour(peptide, index, polar_aa);
  *(features++) = IndexMaxPartialSum(peptide, index, 5);
  *(features++) = IndexMaxPartialSum(peptide, index, 2);
  *(features++) = IndexMinPartialSum(peptide, index, 5);
  *(features++) = IndexMinPartialSum(peptide, index, 2);
  *(features++) = IndexMaxHydrophobicSideHelix(peptide, index);
  *(features++) = IndexMinHydrophobicSideHelix(peptide, index);
  *(features++) = IndexMaxHydrophobicMoment(peptide, index, 100, 11);
  *(features++) = IndexMaxHydrophobicMoment(peptide, index, 180, 11);
  *(features++) = IndexMinHydrophobicMoment(peptide, index, 100, 11);
  *(features++) = IndexMinHydrophobicMoment(peptide, index, 180, 11);
  *(features++) = IndexSumSquaredDiff(peptide, index);
  *(features++) = NumberTypeAA(peptide, polar_aa);
  *(features++) = NumberConsecTypeAA(peptide, polar_aa);
  *(features++) = NumberTypeAA(peptide, hydrophobic_aa);
  *(features++) = NumberConsecTypeAA(peptide, hydrophobic_aa);

  return features;
}

/* get the kPercentageAA*100% AA with the lowest retention and highest retentions*/
pair< set<string>, set<string> > RetentionFeatures::GetExtremeRetentionAA(const map<string, double> &index) {
  vector<string> aa;
  vector<double> retentions;
  map<string, double>::const_iterator it = index.begin();
  int len = index.size();
  int num_aa = (int) ceil(kPercentageAA * len);

  for( ; it != index.end(); ++it) {
    aa.push_back(it->first);
    retentions.push_back(it->second);
  }
  int i;
  bool switched = true;
  double temp_ret;
  string temp_aa;
  while (switched) {
    switched = false;
    for(i = 0; i < len - 1; ++i) {
      if (retentions[i] > retentions[i+1]) {
        switched = true;
        temp_ret = retentions[i];
        retentions[i] = retentions[i+1];
        retentions[i+1] =  temp_ret;
        temp_aa = aa[i];
        aa[i] = aa[i + 1];
        aa[i + 1] = temp_aa;
      }
    }
  }
  // if there are more aa with equal retentions, take all of them
  double val = retentions[num_aa - 1];
  i = num_aa;
  while (i < len && val == retentions[i]) {
    ++i;
  }
  set<string> lowest_set(aa.begin(), aa.begin() + i);
  val = retentions[len - num_aa];
  i = len - num_aa - 1;
  while (i >= 0 && i < len && val == retentions[i]) {
    --i;
  }
  set<string> highest_set(aa.begin() + i + 1, aa.end());

  return pair< set<string>, set<string> >(lowest_set, highest_set);
}

/* calculate the number of a certain type of aa. The set gives the type of these amino acids
 * (for example, the type could be the list "A", "I", "L", "M", "F", "W", "Y", "V", "C" for
 * hydrophobic aa) */
double RetentionFeatures::NumberTypeAA(const string &peptide, const set<string> &amino_acid_type) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  vector<string>::iterator it = amino_acids.begin();
  double occurences = 0.0;

  for( ; it != amino_acids.end(); ++it)
  {
    if (amino_acid_type.find(*it) != amino_acid_type.end()) {
      ++occurences;
    }
  }
  return occurences;
}

/* calculate the number of a consecutibe aa of a certain type. The set gives the type of these such amino acids
 * (for example, the type could be the list "A", "I", "L", "M", "F", "W", "Y", "V", "C" for
 * hydrophobic aa) */
double RetentionFeatures::NumberConsecTypeAA(const string &peptide, const set<string> &amino_acid_type) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  double occurences = 0.0;

  for(int i = 0; i < len - 1; ++i) {
    if ((amino_acid_type.find(amino_acids[i]) != amino_acid_type.end()) &&
        (amino_acid_type.find(amino_acids[i+1]) != amino_acid_type.end()))
      ++occurences;
  }

  return occurences;
}

/* calculate the average hydrophobicity of an index */
double RetentionFeatures::AvgHydrophobicityIndex(const map<string, double> &index) {
  map<string, double>::const_iterator it = index.begin();
  double sum = 0.0;

  for ( ; it != index.end(); ++it)
    sum += it->second;
  return sum / (double) index.size();
}

/* calculate the sum of hydrophobicities of all aa in the peptide; if a modified aa is not in the index
 * then the value of the unmodified aa is used; is an aa is not found at all in the index, the program
 * will abort its execution */
double RetentionFeatures::IndexSum(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  vector<string> amino_acids = GetAminoAcids(peptide);
  vector<string>::iterator it = amino_acids.begin();

  for( ; it != amino_acids.end(); ++it) {
    sum += GetIndexValue(*it, index);
  }
  return sum;
}

/* calculate the average of hydrophobicities of all aa in the peptide; if a modified aa is not in the index
 * then the value of the unmodified aa is used; is an aa is not found at all in the index, the program
 * will abort its execution */
double RetentionFeatures::IndexAvg(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  vector<string> amino_acids = GetAminoAcids(peptide);
  vector<string>::iterator it = amino_acids.begin();

  for( ; it != amino_acids.end(); ++it) {
    sum += GetIndexValue(*it, index);
  }
  return sum / (double) amino_acids.size();
}

/*
 *
 */
/* calculate the hydrophobicity of the N-terminus */
double RetentionFeatures::IndexN(const string &peptide, const map<string, double> &index) {
  int end_position, len = peptide.size();
  string aa;

  aa = peptide.at(0);
  if (aa == "[") {
    end_position = peptide.find("]", 0);
    aa = peptide.substr(0, end_position + 2);
  }
  return GetIndexValue(aa, index);
}

/* calculate the hydrophobicity of the C-terminus */
double RetentionFeatures::IndexC(const string &peptide, const map<string, double> &index) {
  int start_position, len = peptide.size();
  string aa;

  if (peptide.at(len - 2) == ']') {
    start_position = peptide.find_last_of("[", len - 2);
    aa = peptide.substr(start_position);
  } else {
    aa = peptide.at(len - 1);
  }
  return GetIndexValue(aa, index);
}

/* calculate the sum of hydrophobicities of neighbours of polar amino acids */
double RetentionFeatures::IndexNearestNeigbour(const string &peptide, const map<string, double> &index, const set<string> &polar_aa) {
  double sum = 0.0;
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  string aa;

  for(int i = 0; i < len; ++i) {
    aa = amino_acids[i];
    if (polar_aa.find(aa) != polar_aa.end()) {
      if (i > 0) {
        sum += max(0.0, GetIndexValue(amino_acids[i - 1], index));
      }
      if (i < len - 1) {
        sum += max(0.0, GetIndexValue(amino_acids[i + 1], index));
      }
    }
  }
  return sum;
}

/* the most hydrophobic window */
double RetentionFeatures::IndexMaxPartialSum(const string &peptide, const map<string, double> &index, const int &win) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  int window_size = min(win, len - 1);
  double max_sum, sum = 0.0;
  vector<string>::iterator lead = amino_acids.begin();
  vector<string>::iterator lag = amino_acids.begin();

  for( ; lead !=  amino_acids.begin() + window_size; ++lead) {
    sum += GetIndexValue(*lead, index);
  }
  max_sum = sum;
  for( ; lead != amino_acids.end(); ++lead, ++lag) {
    sum -= GetIndexValue(*lag, index);
    sum += GetIndexValue(*lead, index);
    max_sum = max(max_sum, sum);
  }
  return max_sum;
}

/* the least hydrophobic window */
double RetentionFeatures::IndexMinPartialSum(const string &peptide, const map<string, double> &index, const int &win) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  int window_size = min(win, len - 1);
  double min_sum, sum = 0.0;
  vector<string>::iterator lead = amino_acids.begin();
  vector<string>::iterator lag = lead;

  for( ; lead !=  amino_acids.begin() + window_size; ++lead) {
    sum += GetIndexValue(*lead, index);
  }
  min_sum = sum;
  for( ; lead != amino_acids.end(); ++lead, ++lag) {
    sum -= GetIndexValue(*lag, index);
    sum += GetIndexValue(*lead, index);
    min_sum = min(min_sum, sum);
  }
  return min_sum;
}

/* calculate the most hydrophobic sides for alpha helices */
double RetentionFeatures::IndexMaxHydrophobicSideHelix(const string &peptide,  const map<string, double> &index) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  double cos300 = cos(300 * M_PI / 180);
  double cos400 = cos(400 * M_PI / 180);

  if (len < 9) {
    double avg_hydrophobicity_index = AvgHydrophobicityIndex(index);
    return avg_hydrophobicity_index * (1 + 2 * cos300 + 2 * cos400);
  } else {
    double hydrophobicity_side = 0.0;
    double max_hydrophobicity_side = GetIndexValue(amino_acids[4], index) +
                 cos300 * (GetIndexValue(amino_acids[1], index) + GetIndexValue(amino_acids[7], index)) +
                 cos400 * (GetIndexValue(amino_acids[0], index) + GetIndexValue(amino_acids[8], index));

    for(int i = 5; i <= len - 5; ++i) {
      hydrophobicity_side = GetIndexValue(amino_acids[i], index) +
          cos300 * (GetIndexValue(amino_acids[i - 3], index) + GetIndexValue(amino_acids[i + 3], index)) +
          cos400 * (GetIndexValue(amino_acids[i - 4], index) + GetIndexValue(amino_acids[i + 4], index));
      max_hydrophobicity_side = max(max_hydrophobicity_side, hydrophobicity_side);
    }
    return max_hydrophobicity_side;
  }
}

/* calculate the least hydrophobic sides for alpha helices */
double RetentionFeatures::IndexMinHydrophobicSideHelix(const string &peptide,  const map<string, double> &index) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  double cos300 = cos(300 * M_PI / 180);
  double cos400 = cos(400 * M_PI / 180);

  if (len < 9) {
    double avg_hydrophobicity_index = AvgHydrophobicityIndex(index);
    return avg_hydrophobicity_index * (1 + 2 * cos300 + 2 * cos400);
  } else {
    double hydrophobicity_side = 0.0;
    double min_hydrophobicity_side = GetIndexValue(amino_acids[4], index) +
                 cos300 * (GetIndexValue(amino_acids[1], index) + GetIndexValue(amino_acids[7], index)) +
                 cos400 * (GetIndexValue(amino_acids[0], index) + GetIndexValue(amino_acids[8], index));

    for(int i = 5; i <= len - 5; ++i) {
      hydrophobicity_side = GetIndexValue(amino_acids[i], index) +
          cos300 * (GetIndexValue(amino_acids[i - 3], index) + GetIndexValue(amino_acids[i + 3], index)) +
          cos400 * (GetIndexValue(amino_acids[i - 4], index) + GetIndexValue(amino_acids[i + 4], index));
      min_hydrophobicity_side = min(min_hydrophobicity_side, hydrophobicity_side);
    }
    return min_hydrophobicity_side;
  }
}

/* calculate the maximum value of the hydrophobic moment */
double RetentionFeatures::IndexMaxHydrophobicMoment(const string &peptide, const map<string, double> &index, const double &angle_degrees,
                                               const int &win) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  double sin_sum = 0.0, cos_sum = 0.0;
  double angle_radians = angle_degrees * M_PI / 180;

  if (len < win) {
    double avg_hydrophobicity_index = AvgHydrophobicityIndex(index);
    for(int i = 1; i <= win; ++i) {
      cos_sum += cos(i * angle_radians);
      sin_sum += sin(i * angle_radians);
    }
    cos_sum *= avg_hydrophobicity_index;
    sin_sum *= avg_hydrophobicity_index;
    return sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
  } else {
    vector<string>::iterator lead = amino_acids.begin();
    vector<string>::iterator lag = lead;
    double window_hmoment = 0.0;
    int i = 1;
    for( ; lead != amino_acids.begin() + win; ++lead, ++i) {
      cos_sum += GetIndexValue(*lead, index) * cos(i * angle_radians);
      sin_sum += GetIndexValue(*lead, index) * sin(i * angle_radians);
    }
    double max_hmoment = sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
    for( ; lead != amino_acids.end(); ++lead, ++lag, ++i) {
      cos_sum += GetIndexValue(*lead, index) * cos(i * angle_radians);
      cos_sum -= GetIndexValue(*lag, index) * cos((i-win) * angle_radians);
      sin_sum += GetIndexValue(*lead, index) * sin(i * angle_radians);
      sin_sum -= GetIndexValue(*lag, index) * sin((i-win) * angle_radians);
      window_hmoment = sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
      max_hmoment = max(max_hmoment, window_hmoment);
    }
    return max_hmoment;
  }
}

/* calculate the minimum value of the hydrophobic moment */
double RetentionFeatures::IndexMinHydrophobicMoment(const string &peptide, const map<string, double> &index, const double &angle_degrees,
                                               const int &win) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  double sin_sum = 0.0, cos_sum = 0.0;
  double angle_radians = angle_degrees * M_PI / 180;

  if (len < win) {
    double avg_hydrophobicity_index = AvgHydrophobicityIndex(index);
    for(int i = 1; i <= win; ++i) {
      cos_sum += cos(i * angle_radians);
      sin_sum += sin(i * angle_radians);
    }
    cos_sum *= avg_hydrophobicity_index;
    sin_sum *= avg_hydrophobicity_index;
    return sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
  } else {
    vector<string>::iterator lead = amino_acids.begin();
    vector<string>::iterator lag = lead;
    double window_hmoment = 0.0;
    int i = 1;
    for( ; lead != amino_acids.begin() + win; ++lead, ++i) {
      cos_sum += GetIndexValue(*lead, index) * cos(i * angle_radians);
      sin_sum += GetIndexValue(*lead, index) * sin(i * angle_radians);
    }
    double min_hmoment = sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
    for( ; lead != amino_acids.end(); ++lead, ++lag, ++i) {
      cos_sum += GetIndexValue(*lead, index) * cos(i * angle_radians);
      cos_sum -= GetIndexValue(*lag, index) * cos((i-win) * angle_radians);
      sin_sum += GetIndexValue(*lead, index) * sin(i * angle_radians);
      sin_sum -= GetIndexValue(*lag, index) * sin((i-win) * angle_radians);
      window_hmoment = sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum));
      min_hmoment = min(min_hmoment, window_hmoment);
    }
    return min_hmoment;
  }
}

/* calculate the sum of squared differences in hydrophobicities between neighbours */
double RetentionFeatures::IndexSumSquaredDiff(const string &peptide, const map<string, double> &index) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  double squared_diff_sum = 0.0, diff;
  vector<string>::iterator lag = amino_acids.begin();
  vector<string>::iterator lead = lag + 1;

  for( ; lead != amino_acids.end(); ++lead, ++lag) {
    diff = GetIndexValue(*lag, index) - GetIndexValue(*lead, index);
    squared_diff_sum += diff * diff;
  }
  return squared_diff_sum;
}

/* product between hydrophobicity of n- and c- terminus */
/*
double RetentionFeatures::IndexNC(const string &peptide, const map<string, double> &index) {
  double n = max(0.0, IndexN(peptide, index));
  double c = max(0.0, IndexC(peptide, index));
  return n * c;
} */
/* calculate the sum of hydrophobicities of neighbours of D(Aspartic Acid) and E (Glutamic acid) */
/*
double RetentionFeatures::IndexNearestNeigbourNeg(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  string aa;
  char unmodified_aa;

  for(int i = 0; i < len; ++i) {
    aa = amino_acids[i];
    unmodified_aa = GetUnmodifiedAA(aa);
    if (unmodified_aa == 'D' || unmodified_aa == 'E') {
      if (i > 0) {
        sum += max(0.0, GetIndexValue(amino_acids[i - 1], index));
      }
      if (i < len - 1) {
        sum += max(0.0, GetIndexValue(amino_acids[i + 1], index));
      }
    }
  }
  return sum;
} */


/**************************** BULKINESS FUNCTIONS **************************************/
double* RetentionFeatures::ComputeBulkinessFeatures(const string &peptide, const map<string, double> &bulkiness, double *features) {
  *(features++) = ComputeBulkinessSum(peptide, bulkiness);
  return features;
}

/* compute the sum of bulkiness */
double RetentionFeatures::ComputeBulkinessSum(const string &peptide, const map<string, double> &bulkiness) {
  return IndexSum(peptide, bulkiness);
}

/**************************** AMINO ACID FEATURES **************************************/
/* adds a feature giving the number of each of the symbols in the alphabet found in the peptide */
double* RetentionFeatures::FillAAFeatures(const string &peptide, double *retention_features) {
  int number_aa = amino_acids_alphabet_.size();

  for(int i = 0; i < number_aa; ++i) {
    retention_features[i] = 0.0;
  }
  vector<string> amino_acids = GetAminoAcids(peptide);
  vector<string>::iterator it_peptide = amino_acids.begin();
  int i;
  for( ; it_peptide != amino_acids.end(); ++it_peptide) {
    for(i = 0; i < number_aa; ++i) {
      if (amino_acids_alphabet_[i] == (*it_peptide)) {
        ++retention_features[i];
        break;
      }
    }
  }
  return retention_features + amino_acids_alphabet_.size();
}

/**************************** LENGTH FEATURES **************************************/
/* compute the features related to length; */
double* RetentionFeatures::ComputeLengthFeatures(const string &peptide, double *features) {
  *(features++) = PeptideLength(peptide);
  return features;
}

/* calculate the peptide length */
double RetentionFeatures::PeptideLength(const string &peptide) {
  vector<string> amino_acids = GetAminoAcids(peptide);
  return amino_acids.size();
}

/************* FEATURES FOR GROUPS ***************/
/* compute the features when no ptms are present in the data */
double* RetentionFeatures::ComputeNoPTMFeatures(const string &peptide, double *features) {
  // we always compute: kyte and doo little index, svr index, bulkiness, peptide length and amino acid features
  // compute the kyte and doolittle features
  pair< set<string>, set<string> > extreme_aa = GetExtremeRetentionAA(kKyteDoolittle);
  set<string> polar_aa = extreme_aa.first;
  set<string> hydrophobic_aa = extreme_aa.second;
  features = ComputeIndexFeatures(peptide, kKyteDoolittle, polar_aa, hydrophobic_aa, features);

  // compute the svr features
  extreme_aa = GetExtremeRetentionAA(svr_index_);
  polar_aa = extreme_aa.first;
  hydrophobic_aa = extreme_aa.second;
  features = ComputeIndexFeatures(peptide, svr_index_, polar_aa, hydrophobic_aa, features);

  // bulkiness
  features = ComputeBulkinessFeatures(peptide, kBulkiness, features);

  // peptide length
  features = ComputeLengthFeatures(peptide, features);

  // amino acid features
  features = FillAAFeatures(peptide, features);

  return features;
}

// [TO DO: implement this]
/* compute the features when phosphorylations are present in the data */
double* RetentionFeatures::ComputePhosFeatures(const string &peptide, double *features) {
  return features;
}

/************* RETENTION FEATURES FOR PSMS **************/
/* computes the retention features for a set of peptides; return 0 if success */
int RetentionFeatures::ComputeRetentionFeatures(vector<PSMDescription> &psms) {
  vector<PSMDescription>::iterator it= psms.begin();
  for( ; it != psms.end(); ++it) {
    ComputeRetentionFeatures(*it);
  }
  return 0;
}

/* computes the retention features for one psm */
int RetentionFeatures::ComputeRetentionFeatures(PSMDescription &psm) {
  string peptide = psm.getPeptide();
  string::size_type pos1 = peptide.find('.');
  string::size_type pos2 = peptide.find('.', ++pos1);
  string pep = peptide.substr(pos1, pos2 - pos1);
  double* features = psm.getRetentionFeatures();

  // if there is memory allocated
  if (features) {
    if (active_feature_groups_.test(INDEX_NO_PTMS_GROUP)) {
      features = ComputeNoPTMFeatures(pep, features);
    }
    if (active_feature_groups_.test(INDEX_PHOS_GROUP)) {
      features = ComputeNoPTMFeatures(pep, features);
    }
  } else {
    if (VERB >= 1) {
      cerr << "Error: Memory not allocated for the retention features. Execution aborted." << endl;
    }
    exit(1);
  }
  return 0;
}



