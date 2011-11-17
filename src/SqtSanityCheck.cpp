/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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
#include "DataSet.h"
#include "Scores.h"
#include "Globals.h"
#include "SanityCheck.h"
#include "SqtSanityCheck.h"


SqtSanityCheck::SqtSanityCheck(){
}

SqtSanityCheck::~SqtSanityCheck(){
}

const string SqtSanityCheck::fingerPrint = "sqt";

void SqtSanityCheck::getDefaultDirection(vector<vector<double> >& w) {
  // Set init direction to be the most discriminative direction
  int numFeatures = FeatureNames::getNumFeatures();
  for (size_t set = 0; set < w.size(); ++set) {
    for (unsigned int ix = 0; ix < numFeatures + 1; ++ix) {
      w[set][ix] = 0;
    }
    if(numFeatures >= 2) w[set][2] = 1.61; // deltCn
    if(numFeatures >= 3) w[set][3] = 1.1; // Xcorr
    if(numFeatures >= 7) w[set][7] = -0.573; // Peptide length
    if(numFeatures >= 8) w[set][8] = 0.0335; // Charge 1
    if(numFeatures >= 9) w[set][9] = 0.149; // Charge 2
    if(numFeatures >= 10) w[set][10] = -0.156; // Charge 3
  }
  /*

   # first line contains normalized weights, second line the raw weights
   deltCn  Xcorr   PepLen  Charge1 Charge2 Charge3 m0
   1.61    1.1     -0.573  0.0335  0.149   -0.156  -5.23
   20.9    1.99    -0.105  0.328   0.299   -0.312  -9.06

   feat[0]=log(rSp);                     // rank by Sp
   feat[1]=0.0;                     // delt5Cn (leave until last M line)
   feat[2]=0.0;                     // deltCn (leave until next M line)
   feat[3]=xcorr;                   // Xcorr
   feat[4]=sp;                      // Sp
   feat[5]=matched/expected;        // Fraction matched/expected ions
   feat[6]=mass;                    // Observed mass
   feat[7]=peptideLength(pep);      // Peptide length
   feat[8]=(charge==1?1.0:0.0);     // Charge
   feat[9]=(charge==2?1.0:0.0);
   feat[10]=(charge==3?1.0:0.0);
   */
  return;
}

bool SqtSanityCheck::validateDirection(vector<vector<double> >& w) {
  bool ok = SanityCheck::validateDirection(w);
  for (size_t set = 0; set < w.size(); ++set) {
    if (w[set][3] < 0) {
      ok = false;
      if (VERB > 1) {
        cerr << "Warning weight for XCorr negative" << endl;
      }
    }
    if (w[set][2] < 0) {
      ok = false;
      if (VERB > 1) {
        cerr << "Warning weight for deltaCn negative" << endl;
      }
    }
  }
  if (!ok) {
    resetDirection(w);
  }
  return ok;
}

