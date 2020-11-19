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
/* This file defines a number of functions interfacing libsvm */
#include "stdio.h"

#include <vector>
#include <string>

class PSMDescription;
struct svm_parameter;
struct svm_model;

namespace libsvm_wrapper {
  /* train a svr */
  svm_model* TrainModel(const std::vector<PSMDescription*> &psms, const int &number_features, const svm_parameter &parameter);
  /* predict the retention time of psm using the provided svr */
  double PredictRT(const svm_model* svr, const int &number_features, double *features);
  /* save/load a model to/from a file*/
  int SaveModel(FILE* fp, const svm_model* model);
  svm_model* LoadModel(FILE* fp);
}

