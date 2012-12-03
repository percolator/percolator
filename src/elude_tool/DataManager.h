/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@scilifelab.se>

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
#include <utility>

class PSMDescription;

class DataManager {
 public:
   DataManager();
   ~DataManager();
   /* set to null all retention feature pointers and delete the memory allocated for the feature table */
   static void CleanUpTable(std::vector<PSMDescription> &psms, double *feat_table);
   /* load a set of peptides */
   static int LoadPeptides(const std::string &file_name, const bool includes_rt, const bool includes_context,
                           std::vector<PSMDescription> &psms, std::set<std::string> &aa_alphabet);
   /* memory allocation for the feature table; return a pointer to the feature table*/
   static double* InitFeatureTable(const int &no_features, std::vector<PSMDescription> &psms);
   /* remove duplicate peptides */
   static int RemoveDuplicates(std::vector<PSMDescription> &psms);
   /* remove from the train set the peptides that are also in the test set */
   static int RemoveCommonPeptides(const std::vector<PSMDescription> &test_psms,
                            std::vector<PSMDescription> &train_psms);
   /* get the peptide sequence from a peptide given as A.XXX.B */
   static std::string GetMSPeptide(const std::string &peptide);
   /* check if child is in-source fragment of parent*/
   static bool IsFragmentOf(const PSMDescription &child, const PSMDescription &parent,
                     const double &diff, const std::map<std::string, double> &index);
   /* remove in source fragments from the train data; if remove_from_test is true,
    * then we remove the fragments from the test data as well; return a list of
    * in-source fragment, Train/Test depending where the fragment was identified */
   static std::vector< std::pair<PSMDescription, std::string> > RemoveInSourceFragments(
       const double &diff, const std::map<std::string, double> &index,
       bool remove_from_test, std::vector<PSMDescription> &train_psms,
       std::vector<PSMDescription> &test_psms);
   /* combine the train and the test data */
   static std::vector< std::pair<std::pair<PSMDescription, std::string>, bool> > CombineSets(
       std::vector<PSMDescription> &train_psms, std::vector<PSMDescription> &test_psms);
   /* remove non-enzymatic peptides*/
   static std::vector<PSMDescription> RemoveNonEnzymatic(std::vector<PSMDescription> &psms,
       const std::string &mesg);
   /* write in source fragments to file */
   static int WriteInSourceToFile(const std::string &file_name,
       const std::vector< std::pair<PSMDescription, std::string> > &psms);
   /* write peptides to output file */
   static int WriteOutFile(const std::string &file_name,
       const std::vector<PSMDescription> &psms, bool includes_rt);
};

#endif /* DATAMANAGER_H_ */
