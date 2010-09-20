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
 * Class that implement the abstract class SVRModel
 * It provides the required functionalities by using the libsvm package
 */
#ifndef ELUDE_LIBSVRMODEL_H_
#define ELUDE_LIBSVRMODEL_H_

#include "svm.h"
#include "SVRModel.h"

class PSMDescription;

enum SVRType {LINEAR_SVR, RBF_SVR};

class LibSVRModel : SVRModel {
 public:
   LibSVRModel();
   ~LibSVRModel();
   /* calibrate the values of the parameters */
   virtual int CalibrateModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features) {}
   /* train a svr model */
   virtual int TrainModel(const std::vector<PSMDescription> &train_psms, const int &number_features) {}
   /* predict retention time using the trained model */
   virtual int PredictRT(const int &number_features, PSMDescription &psm) {}
   /* save a svr model */
   virtual int SaveModel(const std::ostream &out_stream) {}
   /* load a svr model */
   virtual int LoadModel(const std::istream &input_stream) {}

 private:
   /* grids for parameter calibration; only Epsilon and C are used for linear SVR */
   static const double kGridC[10];
   static const double kGridEpsilon[3];
   static const double kGridGamma[10];
   /* the type of the kernel; could be linear or RBF */
   SVRType kernel_;
   /* svr structure */
   svm_model *svr_;
   /* parameters of the svr */
   svm_parameter svr_parameters_;
};

#endif /* ELUDE_LIBSVRMODEL_H_ */
