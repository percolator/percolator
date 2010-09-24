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
#include <algorithm>

#include "DataManager.h"
#include "RetentionFeatures.h"
#include "PSMDescription.h"
#include "Globals.h"
#include "Enzyme.h"

#define MYABS(x) x >= 0 ? x : (-1) * x

DataManager::DataManager() {
}

DataManager::~DataManager() {
}

/* set to null all retention feature pointers and delete memory */
void DataManager::CleanUpTable(vector<PSMDescription> &psms, double *feat_table) {
  vector<PSMDescription>::iterator it;
  for(it = psms.begin(); it != psms.end(); ++it) {
    it->retentionFeatures = NULL;
  }
  delete[] feat_table;
}

/* load a set of peptides; if the file includes retention time, then includes_rt is true; is the peptides is given in the
 * format A.XXX.B then includes_context is true;  the results are a vector of peptides and a set of all aa present in the peptides */
int DataManager::LoadPeptides(const string &file_name, const bool includes_rt, const bool includes_context,
                              vector<PSMDescription> &psms, set<string> &aa_alphabet) {
  ifstream in(file_name.c_str(), ios::in);
  if (in.fail()) {
    if (VERB >= 1) {
      cerr << "Error: Unable to open " << file_name << ". Execution aborted. " << endl;
    }
    exit(1);
  }
  string peptide_sequence;
  vector<string> amino_acids;
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
double* DataManager::InitFeatureTable(const int &no_features, vector<PSMDescription> &psms) {
  int no_records = psms.size();
  double *feat_pointer = new double[no_records * no_features];

  if (!feat_pointer) {
    if (VERB >= 1) {
      cerr << "Error: Unable to allocate the feature table. Execution aborted."
           << endl;
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

/* remove duplicate peptides */
int DataManager::RemoveDuplicates(std::vector<PSMDescription> &psms) {
  sort(psms.begin(), psms.end(), less<PSMDescription> ());
  psms.resize(distance(psms.begin(), unique(psms.begin(), psms.end())));
  return 0;
}

/* remove from the train set the peptides that are also in the test set */
int DataManager::RemoveCommonPeptides(const vector<PSMDescription> &test_psms,
                                      vector<PSMDescription> &train_psms) {
  vector<PSMDescription>::const_iterator it = test_psms.begin();
  for ( ; it != test_psms.end(); ++it) {
    train_psms.erase(remove(train_psms.begin(), train_psms.end(), (*it)),
                     train_psms.end());
  }
  return 0;
}

// get the peptide information (for A.MCD.E -> return MCD)
string DataManager::GetMSPeptide(const string& peptide) {
  int pos1 = peptide.find('.');
  int pos2 = peptide.find('.', ++pos1);
  return peptide.substr(pos1, pos2 - pos1);
}

/* check if child is in-source fragment from parent either the difference
 * in hydrophobicity is greater than difference according to the givem index*/
bool DataManager::IsFragmentOf(const PSMDescription &child, const PSMDescription &parent,
                               const double &diff, const map<string, double> &index) {
  string peptide_parent = parent.peptide;
  string ms_peptide_parent = GetMSPeptide(peptide_parent);
  string peptide_child = child.peptide;
  string ms_peptide_child = GetMSPeptide(peptide_child);
  // the peptide child has to be included in the larger peptide
  if ((ms_peptide_child.length() >= ms_peptide_parent.length()) ||
      (ms_peptide_parent.find(ms_peptide_child) == string::npos)) {
    return false;
  }
  // parent enzymatic, child non enzymatic
  if (Enzyme::isEnzymatic(peptide_parent) &&
     (!Enzyme::isEnzymatic(peptide_child))) {
    return true;
  }
  // difference in retention sum > diff
  double sum_parent = RetentionFeatures::IndexSum(ms_peptide_parent, index);
  double sum_child = RetentionFeatures::IndexSum(ms_peptide_child, index);
  double retention_difference = sum_child - sum_parent;
  if (MYABS(retention_difference) > diff) {
    return true;
  } else {
    return false;
  }
}

vector< pair<PSMDescription, string> > DataManager::CombineSets(vector<PSMDescription>
       &train_psms, vector<PSMDescription> &test_psms) {
  vector< pair<PSMDescription, string> > temp;
  vector<PSMDescription>::iterator it = train_psms.begin();
  pair<PSMDescription, string> psm_pair;

  for( ; it != train_psms.end(); ++it) {
    psm_pair = make_pair(*it, "train");
    temp.push_back(psm_pair);
  }
  for(it = test_psms.begin() ; it != test_psms.end(); ++it) {
     psm_pair = make_pair(*it, "test");
     temp.push_back(psm_pair);
   }
  return temp;
}

pair<PSMDescription, string> DataManager::RemoveInSourceFragments(
    const double &diff, const map<string, double> &index,
    bool remove_from_test, vector<PSMDescription> &train_psms,
    vector<PSMDescription> &test_psms) {

 vector< pair<PSMDescription, string> > combined_psms;
 combined_psms = CombineSets(train_psms, test_psms);
 // sort the psms according to retention time
 // determine the

}

/*
void RTPredictor::addToPairVector(vector<PSMDescription> psms, bool value,
                                  vector<pair<pair<PSMDescription, bool> ,
                                      bool> > & psmPairs) {
  pair<PSMDescription, bool> tmp1;
  pair<pair<PSMDescription, bool> , bool> tmp2;
  for (int i = 0; i < psms.size(); ++i) {
    tmp1.first = psms[i];
    tmp1.second = value;
    tmp2.first = tmp1;
    tmp2.second = false;
    psmPairs.push_back(tmp2);
  }
}

// remove in source CID peptides
void RTPredictor::removeSourceCIDs(vector<pair<
    pair<PSMDescription, bool> , bool> > & psms) {
  double rt, rtp;
  int noPsms = psms.size(), i, j;
  bool isDecay;
  vector<pair<pair<PSMDescription, bool> , bool> >::iterator it;
  if (VERB > 2) {
    if (removeDecaying) {
      cerr << endl
          << "Removing in source fragments from train and test sets..."
          << endl;
    } else {
      cerr << endl << "Removing in source fragments from train set..."
          << endl;
    }
  }
  sort(psms.begin(), psms.end(), mypair);
  for (i = 0; i < noPsms; ++i) {
    rt = psms[i].first.first.getRetentionTime();
    j = i - 1;
    isDecay = false;
    while ((j >= 0) && (((rtp = psms[j].first.first.getRetentionTime())
        * 1.05) >= rt)) {
      if (isChildOf(psms[i].first.first, psms[j].first.first)) {
        psms[i].second = true;
        isDecay = true;
        break;
      }
      j--;
    }
    if (!isDecay) {
      j = i + 1;
      while ((j < noPsms) && (((rtp
          = psms[j].first.first.getRetentionTime()) * 0.95) <= rt)) {
        if (isChildOf(psms[i].first.first, psms[j].first.first)) {
          psms[i].second = true;
          break;
        }
        j++;
      }
    }
  }
  if (!decayingPeptidesFile.empty()) {
    writeDecayToFile(psms);
  }
  int counts = (int)count_if(psms.begin(), psms.end(), decay);
  if (VERB > 2) {
    cerr << counts << " in source fragments were identified." << endl;
  }
  if (removeDecaying) {
    it = remove_if(psms.begin(), psms.end(), decay);
  } else {
    it = remove_if(psms.begin(), psms.end(), decayTrain);
  }
  psms.resize(distance(psms.begin(), it));
} */
