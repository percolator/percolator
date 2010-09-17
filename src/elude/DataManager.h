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
 * This file stores the class Data Manager, which stores and manipulates the entry data in Elude
 */
#ifndef ELUDE_DATAMANAGER_H_
#define ELUDE_DATAMANAGER_H_

#include <map>
#include <set>
#include <vector>
#include <string>

class PSMDescription;

using namespace std;

class DataManager {
 public:
   DataManager();
   ~DataManager();
   /* load a set of peptides */
   static int LoadPeptides(const string &file_name, const bool includes_rt, const bool includes_context,
                           vector<PSMDescription> &psms, set<string> &aa_alphabet);
   /* memory allocation for the feature table; return a pointer to the feature table*/
   double* InitFeatureTable(const int &no_features, vector<PSMDescription> &psms);

   /************ Accessors and mutators ************/
   inline vector<PSMDescription>& train_psms() { return train_psms_; }
   inline vector<PSMDescription>& test_psms() { return test_psms_; }
   inline set<string>& train_aa_alphabet() { return train_aa_alphabet_; }
   inline set<string>& test_aa_alphabet() { return test_aa_alphabet_; }

 private:
   /* train and test peptide-spectrum matches */
   vector<PSMDescription> train_psms_;
   vector<PSMDescription> test_psms_;
   /* the amino acid alphabet in train and test, respectively */
   set<string> train_aa_alphabet_;
   set<string> test_aa_alphabet_;
   /* pointers to the feature table of the train and test peptides */
   double *train_features_table_, *test_feature_table_;
};

#endif /* DATAMANAGER_H_ */
