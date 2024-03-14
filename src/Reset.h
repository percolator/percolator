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
  Reset();
  ~Reset() { if (pSVMInput_ != nullptr) delete pSVMInput_;};
  int reset(Scores &psms, double selectionFDR, double fractionTrain = 0.5);
  int iterationOfReset(Scores &train, double selectionFDR, double threshold);
  int splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTrain);
 protected:
    AlgIn * pSVMInput_ = nullptr;
};
#endif // RESET_H
