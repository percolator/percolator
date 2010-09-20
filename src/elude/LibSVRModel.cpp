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
 * The file includes definitions of variables and methods in the class LibSVRModel
 */
#include <math.h>

#include "LibSVRModel.h"

/* grids for parameter calibration */
const double LibSVRModel::kGridC[] = {pow(2., -2), pow(2., -1), pow(2., 0), pow(2., 1), pow(2., 2), pow(2., 3), pow(2., 4),
                                pow(2., 5), pow(2., 6), pow(2., 7)};
const double LibSVRModel::kGridEpsilon[] = {0.001, 0.01, 0.1};
const double LibSVRModel::kGridGamma[] = {pow(2., -8), pow(2., -7), pow(2., -6), pow(2., -5), pow(2., -4),
                                    pow(2., -3), pow(2., -2), pow(2., -1), pow(2., 0),  pow(2., 1)};

LibSVRModel::LibSVRModel() : svr_(NULL), kernel_(RBF_SVR) {
}

LibSVRModel::~LibSVRModel() {
}
