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

class PSMDescription;

class DataManager {
 public:
   DataManager();
   ~DataManager();
   /* load a set of peptides; if the file includes retention time, then
    * includes_rt is true; the results are a set of all amino acids present
    * in the peptides and a vector of peptides (last two arguments) */
   static int LoadPeptides(const string &file_name, const bool includes_rt, vector<PSMDescription> &psms, set<string> &aa_alphabet);

   /************ Accessors and mutators ************/
   inline vector<PSMDescription>& train_psms() const { return train_psms_; }
   inline vector<PSMDescription>& test_psms() const { return test_psms_; }
   inline set<string>& train_aa_alphabet() const { return train_aa_alphabet(); }
   inline set<string>& test_aa_alphabet() const { return test_aa_alphabet(); }

 private:
   vector<PSMDescription> train_psms_;
   vector<PSMDescription> test_psms_;
   set<string> train_aa_alphabet;
   set<string> test_aa_alphabet;
};

#endif /* DATAMANAGER_H_ */
