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
 * This file stores the class EludeCaller which is managing the interface with the user
 */

#ifndef ELUDE_ELUDECALLER_H_
#define ELUDE_ELUDECALLER_H_

#include <string>
#include <vector>
#include <set>

#include "PSMDescription.h"

class EludeCaller{
 public:
   EludeCaller();
   ~EludeCaller();
   /* introductory message */
   std::string Greeter() const;
   /* parse the command line arguments */
   bool ParseOptions(int argc, char** argv);

   /************ Accessors and mutators ************/
   inline std::vector<PSMDescription>& train_psms() { return train_psms_; }
   inline std::vector<PSMDescription>& test_psms() { return test_psms_; }
   inline std::set<std::string>& train_aa_alphabet() { return train_aa_alphabet_; }
   inline std::set<std::string>& test_aa_alphabet() { return test_aa_alphabet_; }

 private:
   /* the default library path */
   static std::string library_path_;
   /* lts coverage */
   static double lts_coverage_;
   /* the file including training data */
   std::string train_file_;
   /* the format of the peptides is A.XXX.B; note that here A and B can be '-'
    * The default format is just the peptide sequence */
   bool context_format_;
   /* the file including the test data  */
   std::string test_file_;
   /* Does the test file includes the retention times? By default no */
   bool test_includes_rt_;
   /* file to save the model */
   std::string save_model_file_;
   /* file to load a model from */
   std::string load_model_file_;
   /* the output file */
   std::string output_file_;
   /* select the file from the library? */
   bool automatic_model_sel_;
   /* append the model to the library */
   bool append_model_;
   /* linear calibration? */
   bool linear_calibration_;
   /* remove duplicates from the test set */
   bool remove_duplicates_;
   /* remove in source fragmentation? */
   bool remove_in_source_;
   /* file to save in source fragments */
   std::string in_source_file_;
   /* enzyme */
   std::string the_enzyme_;
   /* remove non enzymatic */
   bool remove_non_enzymatic_;
   /* file to save the retention index */
   std::string index_file_;
   /* remove the peptides from the train set that are also in the test set */
   bool remove_common_peptides_;
   /* train and test peptide-spectrum matches */
   std::vector<PSMDescription> train_psms_;
   std::vector<PSMDescription> test_psms_;
   /* the amino acid alphabet in train and test, respectively */
   std::set<std::string> train_aa_alphabet_;
   std::set<std::string> test_aa_alphabet_;
   /* pointers to the feature table of the train and test peptides */
   double *train_features_table_, *test_features_table_;

};

#endif /* ELUDE_ELUDECALLER_H_ */
