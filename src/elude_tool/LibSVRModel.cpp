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
 * The file includes definitions of variables and methods in the class LibSVRModel
 */
#include <math.h>
#include <stdlib.h>

#include "LibSVRModel.h"
#include "LibsvmWrapper.h"
#include "PSMDescription.h"

using namespace std;

/* grids for parameter calibration */
const double LibSVRModel::kGridC[] = {pow(2., -2), pow(2., -1), pow(2., 0), pow(2., 1), pow(2., 2), pow(2., 3), pow(2., 4),
                                pow(2., 5), pow(2., 6), pow(2., 7)};
const double LibSVRModel::kGridEpsilon[] = {0.001, 0.01, 0.1, 1.0};
const double LibSVRModel::kGridGamma[] = {pow(2., -8), pow(2., -7), pow(2., -6), pow(2., -5), pow(2., -4),
                                    pow(2., -3), pow(2., -2), pow(2., -1), pow(2., 0),  pow(2., 1)};
/* for linear SVR a different grid for C can be used since it gets very slow otherwise */
const double LibSVRModel::kLinearGridC[] = {pow(2., -3), pow(2., -2), pow(2., -1), pow(2., 0), pow(2., 1), pow(2., 2), pow(2., 3)};

/* always 3-fold cross-validation */
const int LibSVRModel::k = 3;

LibSVRModel::LibSVRModel() : svr_(NULL) {
  InitSVRParameters(RBF_SVR);
}

LibSVRModel::LibSVRModel(const SVRType &kernel_type) : svr_(NULL) {
  InitSVRParameters(kernel_type);
}

LibSVRModel::~LibSVRModel() {
  if (svr_) {
    svm_destroy_model(svr_);
    svr_ = NULL;
  }
}

/* initialize the SVR parameter for a RBF kernel*/
int LibSVRModel::InitSVRParameters(const SVRType &kernel_type) {
  svr_parameters_.svm_type = EPSILON_SVR;
  if (kernel_type == LINEAR_SVR) {
    kernel_ = LINEAR_SVR;
    svr_parameters_.kernel_type = LINEAR;
  } else if (kernel_type == RBF_SVR) {
    kernel_ = RBF_SVR;
    svr_parameters_.kernel_type = RBF;
  } else {
    if (VERB >= 4) {
      cerr << "Warning: Unknown kernel type. Using RBF instead. " << endl;
    }
  }
  svr_parameters_.coef0 = 0;
  svr_parameters_.cache_size = 300;
  svr_parameters_.eps = 0.001;
  svr_parameters_.shrinking = 1;
  svr_parameters_.probability = 0;
  svr_parameters_.nr_weight = 0;
  svr_parameters_.weight_label = NULL;
  svr_parameters_.weight = NULL;
  svr_parameters_.C = 0.0;
  svr_parameters_.gamma = 0.0;
  svr_parameters_.p = 0.0;
  return 0;
}

/* set the SVR type */
int LibSVRModel::SetSVRType(const SVRType &type) {
  if (type == LINEAR_SVR) {
    svr_parameters_.kernel_type = LINEAR;
    svr_parameters_.gamma = 0.0;
  } else if (type == RBF_SVR) {
    svr_parameters_.kernel_type = RBF;
  } else {
    if (VERB >= 3) {
      cerr << "Warning: unable to set the kernel type to " << type <<"." << endl;
    }
    return 1;
  }
  return 0;
}

/* set epsilon, C  */
int LibSVRModel::setLinearSVRParam(const double &eps, const double &C){
  svr_parameters_.p = eps;
  svr_parameters_.C = C;
  return 0;
}

/* set epsilon, C, gamma  */
int LibSVRModel::setRBFSVRParam(const double &eps, const double &C, const double &gamma){
  svr_parameters_.p = eps;
  svr_parameters_.C = C;
  svr_parameters_.gamma = gamma;
  return 0;
}

/* train a svr model */
int LibSVRModel::TrainModel(const std::vector<PSMDescription*> &train_psms, const int &number_features) {
  if (svr_) {
    svm_destroy_model(svr_);
  }
  svr_ = libsvm_wrapper::TrainModel(train_psms, number_features, svr_parameters_);
  return 0;
}

/* predict retention time using the trained model */
double LibSVRModel::PredictRT(const int &number_features, double *features) {
  if (svr_) {
    return libsvm_wrapper::PredictRT(svr_, number_features, features);
  }
  else {
    ostringstream temp;
    temp << "Error : No SVR model available. Execution aborted." << endl;
    throw MyException(temp.str());
  }
}

/* predict rt for a set of peptides and return the value of the error */
double LibSVRModel::EstimatePredictionError(const int &number_features, const vector<PSMDescription*> &test_psms) {
  double ms_error = 0.0, predicted_rt = 0.0, deviation;
  vector<PSMDescription*>::const_iterator it = test_psms.begin();

  for ( ; it != test_psms.end(); ++it) {
    predicted_rt =  PredictRT(number_features, (*it)->getRetentionFeatures());
    deviation = predicted_rt - (*it)->getRetentionTime();
    ms_error += deviation * deviation;
  }
  return ms_error / (double)test_psms.size();
}

/* perform k-fold cross validation; return error value */
double LibSVRModel::ComputeKFoldValidation(const std::vector<PSMDescription*> &psms, const int &number_features) {
  vector<PSMDescription*> train, test;
  int len = psms.size();
  // sum of prediction errors
  double sum_pek = 0.0, pek = 0.0;

  for (int i = 0; i < k; ++i) {
    train.clear();
    test.clear();
     // get training and testing sets
     for (int j = 0; j < len; ++j) {
       if ((j % k) == i) {
         test.push_back(psms[j]);
       } else {
         train.push_back(psms[j]);
       }
     }
     TrainModel(train, number_features);
     pek = EstimatePredictionError(number_features, test);
     sum_pek += pek;
   }
  return sum_pek / (double)k;
}

/* calibrate the values of the parameters for a linear SVR; the values of the best parameters
  * are stored in the svr_parameters_ member */
int LibSVRModel::CalibrateLinearModel(const std::vector<PSMDescription*> &calibration_psms,
                                      const int &number_features) {
  double best_c, best_e;
  double best_error = 1e100, error;
  int size_grid_c = sizeof(kLinearGridC) / sizeof(kLinearGridC[0]);
  int size_grid_e = sizeof(kGridEpsilon) / sizeof(kGridEpsilon[0]);
  //int total_steps = grid_c.size() * grid_e.size();
  //int step = 1;

  for(int i = 0; i < size_grid_c; ++i) {
    for(int j = 0; j < size_grid_e; ++j) {
      //++step;
      svr_parameters_.C = kLinearGridC[i];
      svr_parameters_.p = kGridEpsilon[j];
      error = ComputeKFoldValidation(calibration_psms, number_features);
      //cout << "c, epsilon = " << kLinearGridC[i] << ", " << kGridEpsilon[j] << endl;
      //cout << "err = " << error << "\n" << endl;
      if (error < best_error) {
        best_error = error;
        best_c = kLinearGridC[i];
        best_e = kGridEpsilon[j];
      }
    }
  }
  svr_parameters_.C = best_c;
  svr_parameters_.p = best_e;
  return 0;
}

/* calibrate the values of the parameters for a RBF SVR; the values of the best parameters
  * are stored in the svr_parameters_ member */
int LibSVRModel::CalibrateRBFModel(const std::vector<PSMDescription*> &calibration_psms,
                                   const int &number_features) {
  double best_c, best_e, best_g;
  double best_error = 1e100, error;
  int size_grid_c = sizeof(kGridC) / sizeof(kGridC[0]);
  int size_grid_e = sizeof(kGridEpsilon) / sizeof(kGridEpsilon[0]);
  int size_grid_g = sizeof(kGridGamma) / sizeof(kGridGamma[0]);
  //int total_steps = grid_c.size() * grid_e.size();
  //int step = 1;

  for(int i = 0; i < size_grid_c; ++i) {
    for(int j = 0; j < size_grid_e; ++j) {
      for(int p = 0; p < size_grid_g; ++p) {
        //++step;
        svr_parameters_.C = kGridC[i];
        svr_parameters_.p = kGridEpsilon[j];
        svr_parameters_.gamma = kGridGamma[p];
        error = ComputeKFoldValidation(calibration_psms, number_features);
        //cout << "c, epsilon, gamma = " << kGridC[i] << ", " << kGridEpsilon[j] << ", " << kGridGamma[p] << endl;
        //cout << "err = " << error << "\n" << endl;
        if (error < best_error) {
          best_error = error;
          best_c = kGridC[i];
          best_e = kGridEpsilon[j];
          best_g = kGridGamma[p];
        }
      }
    }
  }
  svr_parameters_.C = best_c;
  svr_parameters_.p = best_e;
  svr_parameters_.gamma = best_g;
  //cout << "---------------" << endl;
  //cout << "c, epsilon, gamma = " << best_c << ", " << best_e << ", " << best_g << endl;
  //cout << "err = " << best_error << "\n" << endl;

  return 0;
}

int LibSVRModel::CalibrateModel(const std::vector<PSMDescription*> &calibration_psms,
                   const int &number_features) {
  if (kernel_ == LINEAR_SVR) {
    return CalibrateLinearModel(calibration_psms, number_features);
  } else {
    return CalibrateRBFModel(calibration_psms, number_features);
  }
}

// save a model
int LibSVRModel::SaveModel(FILE *fp) {
  libsvm_wrapper::SaveModel(fp, svr_);
  return 0;
}

// load a model
int LibSVRModel::LoadModel(FILE *fp) {
  if (svr_ != NULL) {
    svm_destroy_model(svr_);
  }
  svr_ = libsvm_wrapper::LoadModel(fp);
  svr_parameters_ = svr_->param;
  int type = svr_parameters_.kernel_type;
  if (type == 0) {
    kernel_ = LINEAR_SVR;
  } else if (type == 2) {
    kernel_ = RBF_SVR;
  } else {
    ostringstream temp;
    temp << "Error: Kernel type " << type << " is not supported. "
           << "Execution aborted." << endl;
    
    throw MyException(temp.str());
  }
}
