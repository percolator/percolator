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
/* This file defines a number of functions interfacing libsvm */
#include <vector>
#include <istream>
#include <ostream>
#include <iostream>
#include "svm.h"

class PSMDescription;

namespace libsvm_wrapper {
  int CalibrateLinearModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features, const double grid_c[],
                           const double grid_epsilon[], svm_parameter *best_parameters) {}
  int CalibrateRBFModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features, const double grid_c[],
                           const double grid_epsilon[], const double grid_gamma[], svm_parameter *best_parameters);
  double ComputeKFoldValidation(const std::vector<PSMDescription> &psms, const int &number_features, const svm_parameter &parameter, const int &k);
  svm_model* TrainModel(const std::vector<PSMDescription> &psms, const int &number_features, const svm_parameter &parameter);
  double PredictRT(const svm_model* svr, const PSMDescription &psm);
  int LoadModel(std::istream & in_stream);
  int SaveModel(std::ostream & out_stream);
}

