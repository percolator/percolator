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
 * This is just an abstract class storing the methods that any svr should provide
 */
#ifndef ELUDE_SVRMODEL_H_
#define ELUDE_SVRMODEL_H_

#include "stdio.h"
#include <vector>
#include <ostream>
#include <istream>

class PSMDescription;

class SVRModel {
 public:
    /* calibrate the values of the parameters */
   virtual int CalibrateModel(const std::vector<PSMDescription> &calibration_psms, const int &number_features) = 0;
   /* train a svr model */
   virtual int TrainModel(const std::vector<PSMDescription> &train_psms, const int &number_features) = 0;
   /* predict retention time using the trained model */
   virtual double PredictRT(const int &number_features, double *features) = 0;
   /* save a svr model */
   virtual int SaveModel(FILE *fp) = 0;
   /* load a svr model */
   virtual int LoadModel(FILE *fp) = 0;
};

#endif /* ELUDE_SVRMODEL_H_ */
