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

class DataManager {
 public:
   DataManager();
   ~DataManager();
   /* set to null all retention feature pointers and delete the memory allocated for the feature table */
   void CleanUpTable(std::vector<PSMDescription> &psms, double *feat_table);
   /* load a set of peptides */
   static int LoadPeptides(const std::string &file_name, const bool includes_rt, const bool includes_context,
                           std::vector<PSMDescription> &psms, std::set<std::string> &aa_alphabet);
   /* memory allocation for the feature table; return a pointer to the feature table*/
   double* InitFeatureTable(const int &no_features, std::vector<PSMDescription> &psms);

   /************ Accessors and mutators ************/
   inline std::vector<PSMDescription>& train_psms() { return train_psms_; }
   inline std::vector<PSMDescription>& test_psms() { return test_psms_; }
   inline std::set<std::string>& train_aa_alphabet() { return train_aa_alphabet_; }
   inline std::set<std::string>& test_aa_alphabet() { return test_aa_alphabet_; }

 private:
   /* train and test peptide-spectrum matches */
   std::vector<PSMDescription> train_psms_;
   std::vector<PSMDescription> test_psms_;
   /* the amino acid alphabet in train and test, respectively */
   std::set<std::string> train_aa_alphabet_;
   std::set<std::string> test_aa_alphabet_;
   /* pointers to the feature table of the train and test peptides */
   double *train_features_table_, *test_features_table_;
};

#endif /* DATAMANAGER_H_ */
