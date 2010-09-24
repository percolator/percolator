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
    RetentionModel();
   ~RetentionModel();
   /* normalize the features of a set of peptides */
   int NormalizeFeatures(std::vector<PSMDescription> &psms);
   /* init a libSVM model; if linear, then a linear kernel is used; else RBF */
   int InitSVR(const bool linear_svr);

   /* Accessors and mutators */
   inline RetentionFeatures& retention_features() { return retention_features_; }
 private:
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
