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

 *******************************************************************************
/*
 * This file implements the methods of the class RTModel
 */
#include "stdio.h"
#include <algorithm>
#include <fstream>

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
int RetentionModel::NormalizeFeatures(const bool set_set,
    std::vector<PSMDescription> &psms) {
  //cout << psms[0] << endl;
  vector<double*> tmp;
  vector<double*> tmp_ret_feat = PSMDescription::getRetFeatures(psms);
  int number_active_features = retention_features_.GetTotalNumberFeatures();
  the_normalizer_->resizeVecs(number_active_features);
  if (set_set) {
    the_normalizer_->setSet(tmp, tmp_ret_feat, 0, number_active_features);
    vsub_ = the_normalizer_->GetVSub();
    vdiv_ = the_normalizer_->GetVDiv();
  } else {
    the_normalizer_->SetSubDiv(vsub_, vdiv_);
  }
  the_normalizer_->normalizeSet(tmp, tmp_ret_feat);
  //cout << psms[0] << endl;
  return 0;
}

void RetentionModel::PrintSub() {
  for(int i = 0; i < vsub_.size(); ++i)
     cout << vsub_[i] << " ";
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
  if (VERB >= 4) {
    cerr << "Training retention index..." << endl;
  }
  if (!normalized_rts) {
    PSMDescription::setPSMSet(psms);
    sub_ = PSMDescription::normSub;
    div_ = PSMDescription::normDiv;
    PSMDescription::normalizeRetentionTimes(psms);
  }
  // set the amino acids alphabet
  vector<string> alphabet(aa_alphabet.begin(), aa_alphabet.end());
  retention_features_.set_amino_acids_alphabet(alphabet);
  // set the aa as the active features
  bitset<RetentionFeatures::NUM_FEATURE_GROUPS> tmp;
  retention_features_.set_active_feature_groups(tmp.set(RetentionFeatures::AA_GROUP));
  // compute the amino acid features
  retention_features_.ComputeRetentionFeatures(psms);
  // normalize the features
  NormalizeFeatures(true, psms);
  // initialize a linear svr
  InitSVR(true);
  // calibrate the model
  int number_features = retention_features_.GetTotalNumberFeatures();
  svr_model_->CalibrateModel(psms, number_features);
  // train the model
  svr_model_->TrainModel(psms, number_features);
  // get the index (weights of the linear SVR)
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
  return GetRetentionIndex();
}

/* train a retention model; the svr_index is given as a parameter */
int RetentionModel::TrainRetentionModel(const set<string> &aa_alphabet, const map<string, double> &index,
    const bool normalized_rts, vector<PSMDescription> &psms) {
  if (VERB >= 4) {
    cerr << "Training retention model..." << endl;
  }
  // normalize
  if (!normalized_rts) {
    PSMDescription::setPSMSet(psms);
    sub_ = PSMDescription::normSub;
    div_ = PSMDescription::normDiv;
    PSMDescription::normalizeRetentionTimes(psms);
  }
  // set the amino acids alphabet and the index
  vector<string> alphabet(aa_alphabet.begin(), aa_alphabet.end());
  retention_features_.set_amino_acids_alphabet(alphabet);
  retention_features_.set_svr_index(index);
  // set the active features; this should be modified as we add more feature groups
  bitset<RetentionFeatures::NUM_FEATURE_GROUPS> tmp;
  if (index.size() <= 20) {
    retention_features_.set_active_feature_groups(
        tmp.set(RetentionFeatures::INDEX_NO_PTMS_GROUP));
  } else {
    retention_features_.set_active_feature_groups(
        tmp.set(RetentionFeatures::INDEX_PHOS_GROUP));
  }
  // compute the amino acid features
  retention_features_.ComputeRetentionFeatures(psms);
  // normalize the features
  NormalizeFeatures(true, psms);
  // save the sub and the div
  vsub_ = the_normalizer_->GetVSub();
  vdiv_ = the_normalizer_->GetVDiv();

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
  // unnormalize the retention time
  PSMDescription::unnormalizeRetentionTimes(psms);
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
  return 0;
}

/* check if the first set is included in the second vector */
bool RetentionModel::IsSetIncluded(const set<string> &alphabet1,
    const vector<string> &alphabet2, const bool ignore_ptms) {
  set<string>::const_iterator it = alphabet1.begin();
  string aa;
  for( ; it!= alphabet1.end(); ++it) {
    if (find(alphabet2.begin(), alphabet2.end(), (*it)) == alphabet2.end()) {
      if (ignore_ptms) {
        if (find(alphabet2.begin(), alphabet2.end(), (*it).substr(0,1))
            == alphabet2.end()) {
          return false;
        }
      } else {
        return false;
      }
    }
  }
  return true;
}

/* return true if the given set is included in the alphabet used in the retention features */
bool RetentionModel::IsIncludedInAlphabet(const set<string> &alphabet,
    const bool ignore_ptms) {
  return  IsSetIncluded(alphabet, retention_features_.amino_acids_alphabet(),
      ignore_ptms);
}

/* predict rt for a vector of psms; it includes calculation of retention features
* but it assumes that the feature table is initialized.
* Note that we cannot predict rt if the model was not trained on same alphabet */
int RetentionModel::PredictRT(const set<string> &aa_alphabet, const bool ignore_ptms,
    const string &text, vector<PSMDescription> &psms) {
  if (VERB >= 4) {
    cerr << "Predicting retention time for the " << text << " ..." << endl;
  }
  // check if the alphabet of the test psms is in included in the alphabet for which
  // the model was built
  //set<string>::const_iterator it1 = aa_alphabet.begin();
  //for( ; it1 != aa_alphabet.end(); ++it1)
  //  cout << *it1 << " ";
  if (!IsSetIncluded(aa_alphabet, retention_features_.amino_acids_alphabet(),
      ignore_ptms)) {
    return 1;
  }
  // if no SVR available, we cannot predict retention time
  if (!svr_model_) {
    if (VERB >= 2) {
      cerr << "Warning: no svr model available to predict retention time. "<< endl;
    }
    return 1;
  }
  // compute the retention features
  retention_features_.set_ignore_ptms(ignore_ptms);
  retention_features_.ComputeRetentionFeatures(psms);
    // normalize the features
  NormalizeFeatures(false, psms);
  double predicted_rt;
  PSMDescription::normDiv = div_;
  PSMDescription::normSub = sub_;
  int number_features = retention_features_.GetTotalNumberFeatures();
  vector<PSMDescription>::iterator it = psms.begin();
  for( ; it != psms.end(); ++it) {
    predicted_rt = svr_model_->PredictRT(number_features, it->retentionFeatures);
    it->predictedTime = PSMDescription::unnormalize(predicted_rt);
  }
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
  return 0;
}

int RetentionModel::SaveModelToFile(const string &file_name) {
  FILE* fp = fopen(file_name.c_str(), "w");
  if (VERB >= 4) {
    cerr << "Saving model to file " << file_name << "..." << endl;
  }
  if (fp == NULL) {
    if (VERB >= 2) {
      cerr << "Warning: Unable to open " << file_name << ". The model "
           <<"will not be saved. " << endl;
    }
    return 1;
  }
  svr_model_->SaveModel(fp);
  // number of features
  int number_features = retention_features_.GetTotalNumberFeatures();
  if (vsub_.size() !=  number_features) {
    if (VERB >= 2) {
      cerr << "Warning: Incorrect model. Number of normalized features is different "
           << "from the number of active features. No model file saved. " << endl;
    }
    fclose(fp);
    remove(file_name.c_str());
    return 1;
  }
  fprintf(fp, "Number_features %d\n", number_features);
  // active features groups
  string tmp = retention_features_.active_feature_groups().
      to_string<char,char_traits<char>,allocator<char> >();
  fprintf(fp, "Active_groups %s\n", tmp.c_str());
  // sub and div to scale retention times
  fprintf(fp, "Sub %lf\n", sub_);
  fprintf(fp, "Div %lf\n", div_);
  // sub and div for the normalizer
  fprintf(fp, "VSub");
  vector<double>::iterator it = vsub_.begin();
  for( ; it != vsub_.end(); ++it) {
    fprintf(fp, " %lf", (*it));
  }
  fprintf(fp, "\n");
  fprintf(fp, "VDiv");
  for(it = vdiv_.begin(); it != vdiv_.end(); ++it) {
    fprintf(fp, " %lf", (*it));
  }
  fprintf(fp, "\n");
  // the svr index
  map<string, double> index = retention_features_.svr_index();
  if (index.size() > 0) {
    map<string, double>::iterator it_map = index.begin();
    fprintf(fp, "Index %d", (int) index.size());
    for( ; it_map != index.end(); ++it_map) {
      fprintf(fp, " %s %lf", it_map->first.c_str(), it_map->second);
    }
    fprintf(fp, "\n");
  }
  // the alphabet
  vector<string> alphabet = retention_features_.amino_acids_alphabet();
  vector<string>::iterator it_vec = alphabet.begin();
  fprintf(fp, "AA_alphabet %d", (int) alphabet.size());
  for( ; it_vec != alphabet.end(); ++it_vec) {
    fprintf(fp, " %s", it_vec->c_str());
  }
  fclose(fp);
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
}

/* load the model from a file */
int RetentionModel::LoadModelFromFile(const std::string &file_name) {
  if (VERB >= 4) {
    cerr << "Loading model from file " << file_name << "..." << endl;
  }
  FILE* fp = fopen(file_name.c_str(), "r");
  if (fp == NULL) {
    ostringstream temp;
    temp << "Error: Unable to open " << file_name << ". Execution aborted" << endl;
    throw MyException(temp.str());
  }
  if (svr_model_) {
    delete svr_model_;
  }
  svr_model_ = new LibSVRModel();
  svr_model_->LoadModel(fp);
  // load the number of features

  char dummy[50];
  int number_features, ret;
  ret = fscanf(fp, "%s %d", dummy, &number_features);
  // active features groups
  char active_groups[RetentionFeatures::NUM_FEATURE_GROUPS];
  ret = fscanf(fp, "%s %s", dummy, active_groups);
  retention_features_.set_active_feature_groups(
      bitset<RetentionFeatures::NUM_FEATURE_GROUPS>(string(active_groups)));
  // load sub and div to scale retention times
  ret = fscanf(fp, "%s %lf", dummy, &sub_);
  ret = fscanf(fp, "%s %lf", dummy, &div_);
  // load vsub and vdiv to scale the features
  ret = fscanf(fp, "%s", dummy);
  vsub_.resize(number_features);
  for(int i = 0; i < number_features; ++i) {
    ret = fscanf(fp, "%lf", &vsub_[i]);
  }
  ret = fscanf(fp, "%s", dummy);
  vdiv_.resize(number_features);
  for(int i = 0; i < number_features; ++i) {
   ret =  fscanf(fp, "%lf", &vdiv_[i]);
  }
  // the retention index
  map<string, double> index;
  int size;
  double val;
  char aa[50];
  ret = fscanf(fp, "%s %d", dummy, &size);
  for(int i = 0; i < size; ++i) {
    ret = fscanf(fp, "%s %lf", aa, &val);
    index[aa] = val;
  }
  retention_features_.set_svr_index(index);
  // the aa alphabet
  vector<string> alphabet;
  ret = fscanf(fp, "%s %d", dummy, &size);
  alphabet.reserve(size);
  for(int i = 0; i < size; ++i) {
    ret = fscanf(fp, "%s", aa);
    alphabet.push_back(aa);
  }
  retention_features_.set_amino_acids_alphabet(alphabet);
  fclose(fp);
  if (number_features != retention_features_.GetTotalNumberFeatures()) {
    ostringstream temp;
    temp << "Error: The number of features used to train the model does not match "
         << "with the current settings. Execution aborted." << endl;
    throw MyException(temp.str());
  }
  
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
}

int RetentionModel::SaveRetentionIndexToFile(const string &file_name) {
  if (VERB >= 4) {
  cerr << "Saving index to " << file_name << "..." << endl;
  }
  ofstream out(file_name.c_str());
  if (out.fail()) {
    if (VERB >= 3) {
      cerr << "Warning: Unable to open " << file_name << ". The retention index "
           << "cannot be stored. " << endl;
    }
    return 1;
  }
  out<< "# File generated by Elude; Retention index given below. " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  map<string, double> index = retention_features_.svr_index();
  map<string, double>::iterator it_map = index.begin();
  for( ; it_map != index.end(); ++it_map) {
    out << it_map->first << " : " << it_map->second << endl;
  }
  out.close();
  if (VERB >= 4) {
    cerr << "Done." << endl << endl;
  }
}
