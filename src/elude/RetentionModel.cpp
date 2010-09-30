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
#include <algorithm>

#include "RetentionModel.h"
#include "PSMDescription.h"
#include "LibSVRModel.h"
#include "Normalizer.h"

RetentionModel::RetentionModel(Normalizer *norm) : the_normalizer_(norm),
    svr_model_(NULL), is_linear_(false) {
}

RetentionModel::~RetentionModel() {
  if (svr_model_) {
    delete svr_model_;
    svr_model_ = NULL;
  }
}

/* normalize the features of a set of peptides */
int RetentionModel::NormalizeFeatures(std::vector<PSMDescription> &psms) {
  int number_active_features = retention_features_.GetTotalNumberFeatures();
  the_normalizer_->resizeVecs(number_active_features);
  vector<double*> tmp;
  vector<double*> tmp_ret_feat = PSMDescription::getRetFeatures(psms);
  the_normalizer_->setSet(tmp, tmp_ret_feat, 0, number_active_features);
  the_normalizer_->normalizeSet(tmp, tmp_ret_feat);
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

/* get the weights of a linear SVR which are equivalent to an index*/
//map<string, double> RetentionModel::GetRetentionIndex() {
map<string, double> RetentionModel::GetRetentionIndex() {
  int no_features = retention_features_.GetTotalNumberFeatures();
  double temp[no_features + 1][no_features];
  for (int i = 0; i < no_features + 1; ++i) {
    for (int j = 0; j < no_features; ++j) {
      if (j == (i - 1)) {
        temp[i][j] = 1.0;
      } else {
        temp[i][j] = 0.0;
      }
    }
  }
  // get weights
  double background = svr_model_->PredictRT(no_features, temp[0]);
  vector<string> aa_alphabet = retention_features_.amino_acids_alphabet();
  map<string, double> svr_index;
  for (int i = 1; i < no_features + 1; ++i) {
    svr_index[aa_alphabet[i-1]] = svr_model_->PredictRT(no_features, temp[i])
        - background;
  }
  return svr_index;
}

/* build the retention index */
map<string, double> RetentionModel::BuildRetentionIndex(const set<string> &aa_alphabet,
    const bool normalized_rts, vector<PSMDescription> &psms) {
  if (!normalized_rts) {
    PSMDescription::setPSMSet(psms);
    sub = PSMDescription::normSub;
    div = PSMDescription::normDiv;
    PSMDescription::normalizeRetentionTimes(psms);
  }
  // set the amino acids alphabet
  vector<string> alphabet(aa_alphabet.begin(), aa_alphabet.end());
  retention_features_.set_amino_acids_alphabet(alphabet);
  // set the aa as the active features
  bitset<NUM_FEATURE_GROUPS> tmp;
  retention_features_.set_active_feature_groups(tmp.set(AA_GROUP));
  // compute the amino acid features
  retention_features_.ComputeRetentionFeatures(psms);
  // normalize the features
  NormalizeFeatures(psms);
  // initialize a linear svr
  InitSVR(true);
  // calibrate the model
  int number_features = retention_features_.GetTotalNumberFeatures();
  svr_model_->CalibrateModel(psms, number_features);
  // train the model
  svr_model_->TrainModel(psms, number_features);
  // get the index (weights of the linear SVR)
  return GetRetentionIndex();
}

/* train a retention model; the svr_index is given as a parameter */
int RetentionModel::TrainRetentionModel(const set<string> &aa_alphabet, const map<string, double> &index,
    const bool normalized_rts, vector<PSMDescription> &psms) {
  // normalize
  if (!normalized_rts) {
    PSMDescription::setPSMSet(psms);
    sub = PSMDescription::normSub;
    div = PSMDescription::normDiv;
    PSMDescription::normalizeRetentionTimes(psms);
  }
  // set the amino acids alphabet and the index
  vector<string> alphabet(aa_alphabet.begin(), aa_alphabet.end());
  retention_features_.set_amino_acids_alphabet(alphabet);
  retention_features_.set_svr_index(index);
  // set the active features; this should be modified as we add more feature groups
  bitset<NUM_FEATURE_GROUPS> tmp;
  if (index.size() <= 20) {
    retention_features_.set_active_feature_groups(tmp.set(INDEX_NO_PTMS_GROUP));
  } else {
    retention_features_.set_active_feature_groups(tmp.set(INDEX_PHOS_GROUP));
  }
  // compute the amino acid features
  retention_features_.ComputeRetentionFeatures(psms);
  // normalize the features
  NormalizeFeatures(psms);
  // initialize a RBF SVR
  if (svr_model_) {
    delete svr_model_;
  }
  InitSVR(false);
  // calibrate the model
  int number_features = retention_features_.GetTotalNumberFeatures();
  //cout << number_features << endl;
  //cout << psms[0] << endl;
  svr_model_->CalibrateModel(psms, number_features);
  // train the model
  svr_model_->TrainModel(psms, number_features);

  return 0;
}

/* check if the first set is included in the second vector */
bool RetentionModel::IsSetIncluded(const set<string> &alphabet1,
    const vector<string> &alphabet2) {
  set<string>::const_iterator it = alphabet1.begin();
  for( ; it!= alphabet1.end(); ++it) {
    if (find(alphabet2.begin(), alphabet2.end(), (*it)) == alphabet2.end()) {
      return false;
    }
  }
  return true;
}

/* predict rt for a vector of psms; it includes calculation of retention features
* but it assumes that the feature table is initialized.
* Note that we cannot predict rt if the model was not trained on same alphabet */
int RetentionModel::PredictRT(const set<string> &aa_alphabet, vector <PSMDescription> &psms) {
  // check if the alphabet of the test psms is in included in the alphabet for which
  // the model was built
  if (!IsSetIncluded(aa_alphabet, retention_features_.amino_acids_alphabet())) {
    if (VERB >= 2) {
      cerr << "Warning: the current model cannot be used to predict retention "
           << "time for the test peptdides. " << endl;
    }
    return 1;
  }
  // if no SVR available, we cannot predict retention time
  if (!svr_model_) {
    if (VERB >= 2) {
      cerr << "Warning: no svr model available to predict retention time. "<< endl;
    }
    return 1;
  }
  // normalize retention times
  retention_features_.ComputeRetentionFeatures(psms);
  // normalize the features
  vector<double*> tmp;
  vector<double*> tmp_ret_feat = PSMDescription::getRetFeatures(psms);
  the_normalizer_->normalizeSet(tmp, tmp_ret_feat);
  double predicted_rt;
  PSMDescription::normDiv = div;
  PSMDescription::normSub = sub;
  int number_features = retention_features_.GetTotalNumberFeatures();
  vector<PSMDescription>::iterator it = psms.begin();
  for( ; it != psms.end(); ++it) {
    predicted_rt = svr_model_->PredictRT(number_features, it->retentionFeatures);
    it->predictedTime = PSMDescription::unnormalize(predicted_rt);
  }
  return 0;
}
