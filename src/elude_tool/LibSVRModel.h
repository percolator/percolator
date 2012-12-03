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
 * Class that implement the abstract class SVRModel
 * It provides the required functionalities by using the libsvm package
 */
#ifndef ELUDE_LIBSVRMODEL_H_
#define ELUDE_LIBSVRMODEL_H_

#include "Globals.h"
#include "svm.h"
#include "SVRModel.h"

class PSMDescription;

class LibSVRModel : public SVRModel {
 public:
   /* grids for parameter calibration; only Epsilon and C are used for linear SVR */
   static const double kGridC[];
   static const double kGridEpsilon[];
   static const double kGridGamma[];
   /* for linear SVR a different grid can be used */
   static const double kLinearGridC[];
   /* k-fold validation; always k = 3 */
   static const int k;
   enum SVRType {LINEAR_SVR = 0, RBF_SVR = 1};
   LibSVRModel();
   LibSVRModel(const SVRType &kernel_type);
   ~LibSVRModel();
   /* initialize the SVR parameter for a RBF kernel*/
   int InitSVRParameters(const SVRType &kernel_type);
   /* set the SVR type */
   int SetSVRType(const SVRType &type);
   /* set epsilon, C */
   int setLinearSVRParam(const double &eps, const double &C);
   /* set epsilon, C, gamma */
   int setRBFSVRParam(const double &eps, const double &C, const double &gamma);
   /* check if the model is null */
   inline bool IsModelNull() const { svr_ == NULL ? true : false; }
   /* train a svr model */
   virtual int TrainModel(const std::vector<PSMDescription> &train_psms,
                          const int &number_features);
   /* predict retention time using the trained model */
   virtual double PredictRT(const int &number_features, double *features);
   /* predict rt for a set of peptides and return the value of the error */
   double EstimatePredictionError(const int &number_features, const std::vector<PSMDescription> &test_psms);
   /* perform k-fold cross validation; return error value */
   double ComputeKFoldValidation(const std::vector<PSMDescription> &psms, const int &number_features);
   /* calibrate the values of the parameters for a linear SVR; the values of the best parameters
    * are stored in the svr_parameters_ member */
   int CalibrateLinearModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features);
   /* calibrate the values of the parameters for a SVR with RBF kernel; the values of the best parameters
    * are stored in the svr_parameters_ member */
   int CalibrateRBFModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features);
   /* calibrate the values of the parameters */
   virtual int CalibrateModel(const std::vector<PSMDescription> &calibration_psms,
                              const int &number_features);
  /* save a svr model */
   virtual int SaveModel(FILE *fp);
   /* load a svr model */
   virtual int LoadModel(FILE *fp);

   /* Accessors and mutators */
   inline svm_parameter svr_parameters() { return svr_parameters_; }

 private:
   /* the type of the kernel; could be linear or RBF */
   SVRType kernel_;
   /* svr structure */
   svm_model *svr_;
   /* parameters of the svr */
   svm_parameter svr_parameters_;
};

#endif /* ELUDE_LIBSVRMODEL_H_ */
