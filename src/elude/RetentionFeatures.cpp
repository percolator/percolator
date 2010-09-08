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
#include <iostream>
#include <algorithm>

#include "boost/assign.hpp"
#include "RetentionFeatures.h"
#include "Globals.h"

using namespace boost::assign;
using namespace std;

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

RetentionFeatures::RetentionFeatures() {
}

RetentionFeatures::~RetentionFeatures() {
}

/* Return the index value of an modified amino acid; if it is not included in the index,
 * then the value of the unmodified one is returned and an warning is returned */
double RetentionFeatures::GetIndexValue(const string &aa, const map<string, double> &index) {
  map<string, double>::const_iterator index_value = index.find(aa);
  string unmodified_aa;

  if (index_value != index.end()) {
    return index_value->second;
  }
  else {
    unmodified_aa = aa[aa.size() - 1];
    index_value = index.find(unmodified_aa);
    if (index_value != index.end()) {
      if (VERB > 4) {
        cerr << "Warning: Could not find the index value for " << aa << ". Use index value for " << unmodified_aa << " instead" << endl;
      }
      return index_value->second;
    }
    else {
      if (VERB > 1) {
        cerr << "Error: Could not find the index value for " << aa << endl;
        cerr << "Execution aborted" << endl;
      }
      exit(1);
    }
  }
}

/* calculate the sum of hydrophobicities of all aa in the peptide; if a modified aa is not in the index
 * then the value of the unmodified aa is used; is an aa is not found at all in the index, the program
 * will abort its execution */
double RetentionFeatures::IndexSum(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  string aa;
  int end_position, i = 0, len = peptide.size();

  while (i < len) {
    aa = peptide[i];
    if (aa == "[") {
      end_position = peptide.find("]", i);
      aa = peptide.substr(i, end_position - i + 2);
    }
    sum += GetIndexValue(aa, index);
    i += aa.size();
  }

  return sum;
}

/* calculate the average of hydrophobicities of all aa in the peptide; if a modified aa is not in the index
 * then the value of the unmodified aa is used; is an aa is not found at all in the index, the program
 * will abort its execution */
double RetentionFeatures::IndexAvg(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  string aa;
  int end_position, i = 0, len = peptide.size();
  double number_aa = 0.0;

  while (i < len) {
    ++number_aa;
    aa = peptide[i];
    if (aa == "[") {
      end_position = peptide.find("]", i);
      aa = peptide.substr(i, end_position - i + 2);
    }
    sum += GetIndexValue(aa, index);
    i += aa.size();
  }

  return sum / number_aa;
}

/* get the amino acids in a peptide */
vector<string> RetentionFeatures::GetAminoAcids(const string &peptide) {
  vector<string> amino_acids;
  int end_position, i = 0, len = peptide.size();
  string aa;

  while (i < len) {
     aa = peptide[i];
     if (aa == "[") {
       end_position = peptide.find("]", i);
       aa = peptide.substr(i, end_position - i + 2);
     }
     amino_acids.push_back(aa);
     i += aa.size();
   }

  return amino_acids;
}

/* calculate the sum of hydrophobicities of neighbours of R(Argenine) and K (Lysine) */
double RetentionFeatures::IndexNearestNeigbourPos(const string &peptide, const map<string, double> &index) {
  double sum = 0.0;
  vector<string> amino_acids = GetAminoAcids(peptide);
  int len = amino_acids.size();
  string aa;
  char unmodified_aa;

  for(int i = 0; i < len; ++i) {
    aa = amino_acids[i];
    unmodified_aa = aa[aa.size() - 1];
    if (unmodified_aa == 'K' || unmodified_aa == 'R') {
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

/* calculate the peptide length */
double RetentionFeatures::PeptideLength(const string &peptide) {
  string aa;
  int end_position, i = 0, len = peptide.size();
  double number_aa = 0.0;

  while (i < len) {
    ++number_aa;
    if (peptide[i] == '[') {
      end_position = peptide.find("]", i);
      i += end_position - i + 2;
    }
    else {
      ++i;
    }
  }

  return number_aa;
}
