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
 * The file includes definitions of variables and methods in the class DataManager
 */
#include <stdio.h>
#include <assert.h>

#include <cstdlib>
#include <fstream>

#include "DataManager.h"
#include "RetentionFeatures.h"
#include "PSMDescription.h"
#include "Globals.h"

DataManager::DataManager() : train_features_table_(NULL), test_features_table_(NULL) {
  std::string basic_aa[] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};
  train_aa_alphabet_.insert(basic_aa, basic_aa + 20);
  test_aa_alphabet_.insert(basic_aa, basic_aa + 20);
}

DataManager::~DataManager() {
  if (train_features_table_) {
    CleanUpTable(train_psms_, train_features_table_);
    train_features_table_ = NULL;
  }
  if (test_features_table_) {
    CleanUpTable(test_psms_, test_features_table_);
    test_features_table_ = NULL;
  }
}

/* set to null all retention feature pointers and delete memory */
void DataManager::CleanUpTable(std::vector<PSMDescription> &psms, double *feat_table) {
  std::vector<PSMDescription>::iterator it;
  for(it = psms.begin(); it != psms.end(); ++it) {
    it->retentionFeatures = NULL;
  }
  delete[] feat_table;
}

/* load a set of peptides; if the file includes retention time, then includes_rt is true; is the peptides is given in the
 * format A.XXX.B then includes_context is true;  the results are a std::vector of peptides and a set of all aa present in the peptides */
int DataManager::LoadPeptides(const std::string &file_name, const bool includes_rt, const bool includes_context,
                              std::vector<PSMDescription> &psms, std::set<std::string> &aa_alphabet) {
  ifstream in(file_name.c_str(), ios::in);
  if (in.fail()) {
    if (VERB >= 1) {
      cerr << "Error: Unable to open " << file_name << ". Execution aborted. " << endl;
    }
    exit(1);
  }
  std::string peptide_sequence;
  std::vector<std::string> amino_acids;
  int len;
  if (includes_rt) {
    double retention_time;
    while (in >> peptide_sequence >> retention_time) {
      psms.push_back(PSMDescription(peptide_sequence, retention_time));
      if (includes_context) {
        len = peptide_sequence.length();
        peptide_sequence = peptide_sequence[0] + peptide_sequence.substr(2, len - 4) + peptide_sequence[len - 1];
      }
      amino_acids = RetentionFeatures::GetAminoAcids(peptide_sequence);
      aa_alphabet.insert(amino_acids.begin(), amino_acids.end());
    }
  } else {
    while (in >> peptide_sequence) {
      psms.push_back(PSMDescription(peptide_sequence, -1.0));
      if (includes_context) {
        len = peptide_sequence.length();
        peptide_sequence = peptide_sequence[0] + peptide_sequence.substr(2, len - 4) + peptide_sequence[len - 1];
      }
      amino_acids = RetentionFeatures::GetAminoAcids(peptide_sequence);
      aa_alphabet.insert(amino_acids.begin(), amino_acids.end());
   }
  }
  in.close();
  return 0;
}

/* memory allocation for the feature table; return a pointer to the feature table*/
double* DataManager::InitFeatureTable(const int &no_features, std::vector<PSMDescription> &psms) {
  int no_records = psms.size();
  double *feat_pointer = new double[no_records * no_features];

  if (!feat_pointer) {
    if (VERB >= 1) {
      cerr << "Error: Unable to allocate memory for the feature table. Execution aborted. " << endl;
    }
    exit(1);
  } else {
    double *ptr = feat_pointer;
    for (int i = 0; i < no_records; i++, ptr += no_features) {
      psms[i].retentionFeatures = ptr;
    }
    return feat_pointer;
  }
}
