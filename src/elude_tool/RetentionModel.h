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
 * This file stores the class RTModel incapsulating data and methods to build and use a
 * retention model
 */

#ifndef ELUDE_RTMODEL_H_
#define ELUDE_RTMODEL_H_

#include <vector>

#include "RetentionFeatures.h"
#include "SVRModel.h"

class Normalizer;
class PSMDescription;

class RetentionModel {
 public:
    RetentionModel(Normalizer *normalizer);
   ~RetentionModel();
   /* normalize the features of a set of peptides */
   int NormalizeFeatures(const bool set_set, std::vector<PSMDescription*> &psms);
   /* init a libSVM model; if linear, then a linear kernel is used; else RBF */
   int InitSVR(const bool linear_svr);
   /* set the alphabet */
   inline void SetAlphabet(const std::vector<std::string> alphabet) {
       retention_features_.set_amino_acids_alphabet(alphabet); }
   /* build a retention index from a set of psms
    * assume the rt of the psms is already normalized */
   std::map<std::string, double>  BuildRetentionIndex(
       const std::set<std::string> &aa_alphabet, const bool normalized_rts,
       std::vector<PSMDescription*> &psms);
   /*  get the weights of a linear SVR */
   std::map<std::string, double> GetRetentionIndex();
   /* train a retention model; the svr_index is given as a parameter */
   int TrainRetentionModel(const std::set<std::string> &aa_alphabet,
       const std::map<std::string, double> &index, const bool normalized_rts,
       std::vector<PSMDescription*> &psms);
   /* return true if the model is null */
   bool IsModelNull() const { return svr_model_ == NULL; }
   /* check if the first set is included in the second vector */
   static bool IsSetIncluded(const std::set<std::string> &alphabet1,
       const std::vector<std::string> &alphabet2, const bool ignore_ptms);
   /* predict rt for a vector of psms; it includes calculation of retention features
    * but it assumes that the feature table is initialized */
   int PredictRT(const std::set<std::string> &aa_alphabet, const bool ignore_ptms,
       const std::string &text, std::vector<PSMDescription*> &psms);
   /* save the model to a file */
   int SaveModelToFile(const std::string &file_name);
   /* load the model from a file */
   int LoadModelFromFile(const std::string &file_name);
   /* save the retention index to a file */
   int SaveRetentionIndexToFile(const std::string &file_name);
   /* print vsub_ */
   void PrintSub();
   /* return true if the given set is included in the alphabet used in
    * the retention features */
   bool IsIncludedInAlphabet(const std::set<std::string> &alphabet,
       const bool ignore_ptms);

   /* Accessors and mutators */
   inline RetentionFeatures& retention_features() { return retention_features_; }
   inline double sub() const { return sub_; }
   inline double div() const { return div_; }
   inline void set_sub(const double &s) { sub_ = s; }
   inline void set_div(const double &d) { div_ = d; }
   inline void set_index(const std::map<std::string, double> &index) {
	   retention_features_.set_svr_index(index); }

 private:
   /* used to normalize retention times and features*/
   double sub_, div_;
   std::vector<double> vsub_, vdiv_;
   /* the SVR model */
   SVRModel *svr_model_;
   /* is the SVR linear? */
   bool is_linear_;
   /* the normalizer */
   Normalizer *the_normalizer_;
   /* the retention features */
   RetentionFeatures retention_features_;
};

#endif /* ELUDE_RTMODEL_H_ */
