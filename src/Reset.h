/*******************************************************************************
 Copyright 2024 Lukas Käll <lukas.kall@scilifelab.se>

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
 
#ifndef RESET_H
#define RESET_H

#include "ssl.h"
class Scores;

class Reset {
 public:
    Reset() {
        options_.lambda = 1.0;
        options_.lambda_u = 1.0;
        options_.epsilon = EPSILON;
        options_.cgitermax = CGITERMAX;
        options_.mfnitermax = MFNITERMAX;
    };
    ~Reset() { if (pSVMInput_ != nullptr) delete pSVMInput_;};
    int reset(Scores &psms, double selectionFDR, double fractionTraining = 0.5, int decoysPerTarget = 1);
    int iterationOfReset(Scores &train, double selectionFDR, double threshold);
    int splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTrain);

 protected:
    AlgIn * pSVMInput_ = nullptr;
    std::vector<double> w_; // linear scoring weights from SVM
    options options_;

};
#endif // RESET_H