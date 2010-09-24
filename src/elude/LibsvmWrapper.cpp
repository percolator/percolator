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
/* This file includes the implementations for the functions defined in LibsvmWrapper.h */
#include <stdlib.h>

#include "LibsvmWrapper.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "svm.h"

svm_model* libsvm_wrapper::TrainModel(const std::vector<PSMDescription> &psms, const int &number_features, const svm_parameter &parameter) {
  svm_model *svr_model;
  int number_examples = psms.size();
  svm_problem data;
  data.l = number_examples;
  data.x = new svm_node[number_examples];
  data.y = new double[number_examples];
  for (int i = 0; i < number_examples; i++) {
    data.x[i].values = psms[i].retentionFeatures;
    data.x[i].dim = number_features;
    data.y[i] = psms[i].retentionTime;
  }
  // build a model by training the SVM on the given training set
  char const *error_message = svm_check_parameter(&data, &parameter);
  if (error_message != NULL) {
    if (VERB >= 1) {
      cerr << "ERROR: Incorrect parameters for the SVR. Execution aborted. " << endl;
    }
    exit(1);
  }
  svr_model = svm_train(&data, &parameter);
  delete[] data.x;
  delete[] data.y;
  return svr_model;
}

double libsvm_wrapper::PredictRT(const svm_model* svr, const int &number_features, double *features) {
  svm_node node;
  node.values = features;
  node.dim = number_features;

  return svm_predict(svr, &node);
}
