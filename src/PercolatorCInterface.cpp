/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"
#include "PercolatorCInterface.h"

static Caller* pCaller = NULL;
static NSet nset = FOUR_SETS;
static unsigned int numFeatures = 0;
static SetHandler::Iterator* normal = NULL;
static SetHandler::Iterator* decoy1 = NULL;
static SetHandler::Iterator* decoy2 = NULL;
static SetHandler::Iterator* decoy3 = NULL;

Caller* getCaller() {
  if (pCaller == NULL) {
    cerr << "Object pCaller not properly assigned" << endl;
    exit(-1);
  }
  return pCaller;
}

/** Call that initiates percolator */
void pcInitiate(NSet sets, unsigned int numFeat, unsigned int numSpectra,
                char** featureNames, double pi0) {
  pCaller = new Caller();
  nset = sets;
  numFeatures = numFeat;
  pCaller->filelessSetup(numFeatures, numSpectra, featureNames, pi0);
  normal = new SetHandler::Iterator(pCaller->getSetHandler(), Caller::NORMAL);
  decoy1 = new SetHandler::Iterator(pCaller->getSetHandler(), Caller::SHUFFLED);
  if (nset > 2) {
    cerr
      << "This version of percolator only supports 1 decoy set. Percolator was called with nset="
      << nset << endl;
  }
}

/** Call that sets verbosity level
 *  0 is quiet, 2 is default, 6 is more than you want */
void pcSetVerbosity(int verbosity) {
  if (verbosity < 0 || verbosity > 6) {
    cerr
        << "Wrong value of verbosity, should be between 0 and 6, you tried "
        << verbosity << endl;
    exit(-1);
  }
  Globals::getInstance()->setVerbose(verbosity);
}

/** Register a PSM */
void pcRegisterPSM(SetType set, char* identifier, double* features) {
  if ((int)set > (int)nset) {
    cerr << "Tried to access undefined set" << endl;
    exit(-1);
  }
  double* vec = NULL;
  switch (set) {
    case TARGET:
      vec = normal->getNext()->features;
      break;
    case DECOY1:
      vec = decoy1->getNext()->features;
      break;
    case DECOY2:
      vec = decoy2->getNext()->features;
      break;
    case DECOY3:
      vec = decoy3->getNext()->features;
      break;
  }
  if (vec == NULL) {
    cerr << "Pointer out of bound" << endl;
    exit(-1);
  }
  for (unsigned int ix = 0; ix < numFeatures; ix++) {
    vec[ix] = features[ix];
  }
}

/** Function called when we want to start processing */
void pcExecute() {
  pCaller->fillFeatureSets();
  vector<vector<double> > w(1,
                            vector<double> (FeatureNames::getNumFeatures()
                                + 1));
  pCaller->preIterationSetup(w);
  pCaller->train(w);
  pCaller->getFullSet()->calcScores(w[0]);
}

/** Function called when retrieving target scores and q-values after processing,
 * the array should be numSpectra long and will be filled in the same order
 * as the features were inserted */
void pcGetScores(double* scoreArr, double* qArr) {
  int ix = 0;
  SetHandler::Iterator iter(pCaller->getSetHandler(), Caller::NORMAL);
  while (PSMDescription* pPSM = iter.getNext()) {
    ScoreHolder* pSh =
        pCaller->getFullSet()->getScoreHolder(pPSM->features);
    scoreArr[ix] = pSh->score;
    qArr[ix++] = pPSM->q;
  }
}

/** Function that should be called after processing finished */
void pcCleanUp() {
  if (pCaller) {
    delete pCaller;
    pCaller = NULL;
  }
  if (normal) {
    delete normal;
    normal = NULL;
  }
  if (decoy1) {
    delete decoy1;
    decoy1 = NULL;
  }
  if (decoy2) {
    delete decoy2;
    decoy2 = NULL;
  }
  if (decoy3) {
    delete decoy3;
    decoy3 = NULL;
  }
  Globals::clean();
}
