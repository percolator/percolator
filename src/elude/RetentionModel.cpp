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
 * This file implements the methods of the class RTModel
 */
#include "RetentionModel.h"
#include "PSMDescription.h"
#include "LibSVRModel.h"
#include "Normalizer.h"

RetentionModel::RetentionModel() : the_normalizer_(NULL), svr_model_(NULL), is_linear_(false) {
}

RetentionModel::~RetentionModel() {
  if (the_normalizer_) {
    delete the_normalizer_;
    the_normalizer_ = NULL;
  }
  if (svr_model_) {
    delete svr_model_;
    svr_model_ = NULL;
  }
}

/* normalize the features of a set of peptides */
int RetentionModel::NormalizeFeatures(std::vector<PSMDescription> &psms) {
  int number_active_features = retention_features_.GetTotalNumberFeatures();
  the_normalizer_ = Normalizer::getNormalizer();
  the_normalizer_->resizeVecs(number_active_features);
  vector<double*> tmp;
  vector<double*> tRetFeat = PSMDescription::getRetFeatures(psms);
  the_normalizer_->setSet(tmp, tRetFeat, 0, number_active_features);
  the_normalizer_->normalizeSet(tmp, tRetFeat);
  return 0;
}

/* init a libSVM model; if linear, then a linear kernel is used; else RBF */
int RetentionModel::InitSVR(const bool linear_svr) {
  if (linear_svr) {
    svr_model_ = new LibSVRModel(LibSVRModel::LINEAR_SVR);
    is_linear_ = true;
  } else {
    svr_model_ = new LibSVRModel(LibSVRModel::RBF_SVR);
    is_linear_ = false;
  }
}
