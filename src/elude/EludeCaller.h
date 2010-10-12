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
#include "RetentionModel.h"

class Normalizer;

class EludeCaller{
 public:
   EludeCaller();
   ~EludeCaller();
   /* fraction of the peptides wwhen computing the window */
   static const double kFractionPeptides;
   /* introductory message */
   std::string Greeter() const;
   /* parse the command line arguments */
   bool ParseOptions(int argc, char** argv);
   /* set the enzyme type */
   void SetEnzyme(const string &enzyme);
   /* process the train data */
   int ProcessTrainData();
   /* normalize retention times for a set of peptides */
   int NormalizeRetentionTimes(vector<PSMDescription> &psms);
   /* train a retention model */
   int TrainRetentionModel();
   /* process the test data */
   int ProcessTestData();
   /* Load the best model from the library; the function returns a pair consisting of
    * the index of this model in the vector of models and the rank correlation
    * obtained on the calibration peptides using this model */
   std::pair<int, double> AutomaticModelSelection();
   /* main function of Elude */
   int Run();
   /*compute Delta t(95%) window */
   static double ComputeWindow(vector<PSMDescription> &psms);
   /* calculate Spearman's rank correlation */
   static double ComputeRankCorrelation(vector<PSMDescription> &psms);
   /* Compute Pearson's correlation coefficient */
   static double ComputePearsonCorrelation(vector<PSMDescription> & psms);
   /* Delete all the RT model */
   void DeleteRTModels();
   /* Return a list of files in a directory */
   static std::vector<std::string> ListDirFiles(const std::string &dir_name);
   /* Check if a file exists */
   static bool FileExists(const string &file);

   /************ Accessors and mutators ************/
   inline std::vector<PSMDescription>& train_psms() { return train_psms_; }
   inline std::vector<PSMDescription>& test_psms() { return test_psms_; }
   inline std::set<std::string>& train_aa_alphabet() { return train_aa_alphabet_; }
   inline std::set<std::string>& test_aa_alphabet() { return test_aa_alphabet_; }
   inline void set_save_model_file(const string &file) { save_model_file_ = file; }
   inline void set_load_model_file(const string &file) { load_model_file_ = file; }
   inline void set_train_file(const string &file) { train_file_ = file; }
   inline void set_test_file(const string &file) { test_file_ = file; }
   inline void set_in_source_file(const string &file) { in_source_file_ = file; }
   inline void set_context_format(bool cf) { context_format_ = cf; }
   inline void set_test_includes_rt(bool tir) { test_includes_rt_ = tir; }
   inline void set_remove_duplicates(bool rd) { remove_duplicates_ = rd; }
   inline void set_remove_in_source(bool rin) { remove_in_source_ = rin; }
   inline void set_non_enzymatic(bool rnz) { remove_non_enzymatic_ = rnz; }
   inline void set_remove_common_peptides(bool rcp) { remove_common_peptides_ = rcp; }
   inline static void set_lib_path(const string &path) { library_path_ = path; }
   inline void set_automatic_model_sel(const bool ams) { automatic_model_sel_ = ams; }

 private:
   /* the default library path */
   static std::string library_path_;
   /* lts coverage */
   static double lts_coverage_;
   /* difference in hydrophobicity to be called in source fragments */
   static double hydrophobicity_diff_;
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
   /* remove in source fragmentation from test? */
   bool remove_in_source_;
   /* file to save in source fragments */
   std::string in_source_file_;
   /* remove non enzymatic */
   bool remove_non_enzymatic_;
   /* file to save the retention index */
   std::string index_file_;
   /* remove the peptides from the train set that are also in the test set */
   bool remove_common_peptides_;
   /* ignore the ptms */
   bool ignore_ptms_;
   /* train and test peptide-spectrum matches */
   std::vector<PSMDescription> train_psms_;
   std::vector<PSMDescription> test_psms_;
   /* the amino acid alphabet in train and test, respectively */
   std::set<std::string> train_aa_alphabet_;
   std::set<std::string> test_aa_alphabet_;
   /* pointers to the feature table of the train and test peptides */
   double *train_features_table_, *test_features_table_;
   /* true if the test data was processed */
   bool processed_test_;
   /* the retention models */
   vector<RetentionModel*> rt_models_;
   /* the retention model */
   RetentionModel *rt_model_;
   std::map<std::string, double> retention_index_;
   /* the normalizer */
   Normalizer *the_normalizer_;
};

#endif /* ELUDE_ELUDECALLER_H_ */
