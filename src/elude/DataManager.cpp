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

/* Class that stores predicates used in different sorting/partition functions */
class Utilities {
public:
 static bool ComparePairs(const pair<pair<PSMDescription, string>, bool> &psm1,
                  const pair<pair<PSMDescription, string>, bool> &psm2) {
  return psm1.first.first.getRetentionTime() < psm2.first.first.getRetentionTime();
}
 static bool IsInSource(const pair<pair<PSMDescription, string>, bool> &psm) {
   return psm.second;
 }

 static bool IsInSourceAndTrain(const pair<pair<PSMDescription, string>, bool> &psm) {
   return psm.second && psm.first.second == "train";
 }
 static bool IsInTrain(const pair<pair<PSMDescription, string>, bool> &psm) {
   return psm.first.second == "train";
 }

 static PSMDescription GetPSM(const pair<pair<PSMDescription, string>, bool> &psm) {
   return psm.first.first;
 }

 static pair<PSMDescription, string> GetPair(const pair<pair<PSMDescription, string>, bool> &psm) {
   return psm.first;
 }

 static bool IsEnzymatic(const PSMDescription &psm) {
   return Enzyme::isEnzymatic(psm.peptide);
 }
};

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
  if (VERB >= 4) {
    cerr << "Loading file " << file_name << "..." << endl;
  }
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
  if (VERB >= 4) {
    cerr << psms.size() << " peptides loaded." << endl << endl;
  }
  return 0;
}

/* memory allocation for the feature table; return a pointer to the feature table*/
double* DataManager::InitFeatureTable(const int &no_features, vector<PSMDescription> &psms) {
  int no_records = psms.size();
  if (VERB >= 4) {
    cerr << "Initializing feature table for " << no_records << "records and "
         << no_features << " features..." << endl;
  }
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
    if (VERB >= 4) {
      cerr << "Done." << endl << endl;
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
  if (VERB >= 4) {
    cerr << "Removing from the train set the peptides featuring in the test set..." << endl;
  }
  int initial_number = train_psms.size();
  if (train_psms.empty()) {
    if (VERB >= 4) {
      cerr << "Warning: no train data available, thus no common peptides "
           <<"between train and test sets." << endl;
    }
    return 0;
  }
  vector<PSMDescription>::const_iterator it = test_psms.begin();
  for ( ; it != test_psms.end(); ++it) {
    train_psms.erase(remove(train_psms.begin(), train_psms.end(), (*it)),
                     train_psms.end());
  }
  if (VERB >= 4) {
    cerr << (train_psms.size() - initial_number) << " peptides were removed."
         << endl << endl;
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
  if (peptide_parent != ms_peptide_parent && peptide_child != ms_peptide_child &&
      Enzyme::isEnzymatic(peptide_parent) && (!Enzyme::isEnzymatic(peptide_child))) {
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

vector< pair<pair<PSMDescription, string>, bool> > DataManager::CombineSets(vector<PSMDescription>
       &train_psms, vector<PSMDescription> &test_psms) {
  vector< pair<pair<PSMDescription, string> , bool> > temp;
  vector<PSMDescription>::iterator it = train_psms.begin();
  pair<PSMDescription, string> psm_pair;

  for( ; it != train_psms.end(); ++it) {
    psm_pair = make_pair(*it, "train");
    temp.push_back(make_pair(psm_pair, false));
  }
  for(it = test_psms.begin() ; it != test_psms.end(); ++it) {
     psm_pair = make_pair(*it, "test");
     temp.push_back(make_pair(psm_pair, false));
   }
  return temp;
}

vector< pair<PSMDescription, string> > DataManager::RemoveInSourceFragments(
  const double &diff, const map<string, double> &index,
  bool remove_from_test, vector<PSMDescription> &train_psms,
  vector<PSMDescription> &test_psms) {
  if (VERB >= 4) {
    if (remove_from_test) {
      cerr << "Removing in-source fragments from test and train data..." << endl;
    } else {
      cerr << "Removing in-source fragments from the train data..." << endl;
    }
  }
  // store information about whether the peptides is in train or test
  vector< pair<pair<PSMDescription, string>, bool> > combined_psms =
      CombineSets(train_psms, test_psms);
  // sort the psms according to retention time
  sort(combined_psms.begin(), combined_psms.end(), Utilities::ComparePairs);
  int number_psms = combined_psms.size();

  // check in source fragmentation
  double rt_child, rt_parent;
  bool is_in_source;
  int i, j;
  for (i = 0; i < number_psms; ++i) {
    rt_child = combined_psms[i].first.first.getRetentionTime();
    j = i - 1;
    is_in_source = false;
    while ((j >= 0) && (((rt_parent =
        combined_psms[j].first.first.getRetentionTime()) * 1.05) >= rt_child)) {
      if (IsFragmentOf(combined_psms[i].first.first, combined_psms[j].first.first, diff, index)) {
        combined_psms[i].second = true;
        is_in_source = true;
        break;
      }
      j--;
    }
    if (!is_in_source) {
      j = i + 1;
      while ((j < number_psms) && (((rt_parent =
          combined_psms[j].first.first.getRetentionTime()) * 0.95) <= rt_child)) {
        if (IsFragmentOf(combined_psms[i].first.first, combined_psms[j].first.first, diff, index)) {
          combined_psms[i].second = true;
          break;
        }
        j++;
      }
    }
  }
  // partition the PSMs according to whether they are in source fragments or not
  vector< pair<pair<PSMDescription, string>, bool> >::iterator it1 =
      partition(combined_psms.begin(), combined_psms.end(), Utilities::IsInSource);
  // save the fragments
  vector< pair<PSMDescription, string> > fragments;
  vector< pair<pair<PSMDescription, string>, bool> >::iterator it2 =
      combined_psms.begin();
  for( ; it2 != it1; ++it2) {
    fragments.push_back(it2->first);
  }
  if (!remove_from_test) {
    it1 = partition(combined_psms.begin(), it1, Utilities::IsInSourceAndTrain);
  }
  // remove in source fragments
  combined_psms.erase(combined_psms.begin(), it1);
  // remove fragments from initial sets
  it1 = partition(combined_psms.begin(), combined_psms.end(), Utilities::IsInTrain);
  transform(combined_psms.begin(), it1, train_psms.begin(), Utilities::GetPSM);
  train_psms.resize(distance(combined_psms.begin(), it1));
  transform(it1, combined_psms.end(), test_psms.begin(), Utilities::GetPSM);
  test_psms.resize(distance(it1, combined_psms.end()));
  if (VERB >= 4) {
    cerr << fragments.size() << " in-source fragments were identified." << endl;
    if (remove_from_test) {
      cerr << "The train set includes now " << train_psms.size() << " peptides." << endl;
      cerr << "The test set includes now " << test_psms.size() << " peptides." << endl << endl;
    } else {
      cerr << "The train set includes now " << train_psms.size() << " peptides." << endl << endl;
    }
  }
  return fragments;
}

/* return a list of non-enzymatic peptides; this peptides are removed from the psms */
vector<PSMDescription> DataManager::RemoveNonEnzymatic(vector<PSMDescription> &psms,
    const string &mesg) {
  if (VERB >= 4) {
    cerr << "Removing non enzymatic peptides from the " << mesg << "..." << endl;
  }
  int initial_size = psms.size();
  vector<PSMDescription>::iterator it =
      partition(psms.begin(), psms.end(), Utilities::IsEnzymatic);
  vector<PSMDescription> non_enzymatic(it, psms.end());
  //for(int i = 0; i < non_enzymatic.size(); ++i)
  //   cout << non_enzymatic[i].peptide << endl;
  psms.resize(distance(psms.begin(), it));
  if (VERB >= 4) {
    cerr << non_enzymatic.size() << " non-enzymatic peptides were identified " << endl;
    cerr << "The " << mesg << " includes now " << psms.size() << " peptides." << endl << endl;
  }
  return non_enzymatic;
}

/* write a list of in source fragmentation to a file */
int DataManager::WriteInSourceToFile(const string &file_name,
    const vector< pair<PSMDescription, string> > &psms) {
  if (VERB >= 4) {
    cerr << "Writing in-source fragments to " << file_name << "..." << endl;
  }
  ofstream out(file_name.c_str());
  if (out.fail()) {
    if (VERB >= 2) {
      cerr << "Warning: Unable to open " << file_name << ". In-source fragments "
           << "cannot be stored. " << endl;
    }
    return 1;
  }
  out<< "# File generated by Elude; In source fragments listed below. " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  out<< "Peptide\tobserved_retention_time\tSet" << endl;
  vector< pair<PSMDescription, string> >::const_iterator it = psms.begin();
  for ( ; it != psms.end(); ++it) {
    out << it->first.peptide << "\t" << it->first.retentionTime << "\t"
        << it->second << endl;
  }
  out.close();
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
  return 0;
}

/* write a set of peptides to an output file */
int DataManager::WriteOutFile(const string &file_name,
    const vector<PSMDescription> &psms, bool includes_rt) {
  if (VERB >= 4) {
    cerr << "Writing predictions to file " << file_name << "..." << endl;
  }
  ofstream out;
  out.open(file_name.c_str());
  if (out.fail()) {
    if (VERB >= 2) {
      cerr << "Warning: Unable to open " << file_name << ". The output file cannot "
           << "be generated. " << endl;
    }
    return 1;
  }
  out << "# File generated by Elude. Predicted retention times are given below. " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  if (includes_rt) {
    out<< "Peptide\tPredicted_RT\tObserved_RT" << endl;
  } else {
    out<< "Peptide\tPredicted_RT" << endl;
  }
  vector<PSMDescription>::const_iterator it = psms.begin();
  for ( ; it != psms.end(); ++it) {
    if (includes_rt) {
      out << it->peptide << "\t" << it->predictedTime << "\t"
          << it->retentionTime << endl;
    } else {
      out << it->peptide << "\t" << it->predictedTime << endl;
    }
  }
  out.close();
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
  return 0;
}
