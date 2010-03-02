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
#include <assert.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "UniNormalizer.h"
#include "Globals.h"

Normalizer::Normalizer() {
}

Normalizer::~Normalizer() {
}

void Normalizer::normalizeSet(set<DataSet *> & setVec) {
  double * features;
  PSMDescription* pPSM;
  set<DataSet *>::iterator it;
  if (VERB > 4) {
    cerr << "First 10 feature vectors before normalization" << endl;
    cerr.precision(3);
    it = setVec.begin();
    int ixPos = -1;
    cerr << "Label of this set is " << (*it)->getLabel() << endl;
    while ((pPSM = (*it)->getNext(ixPos)) != NULL && ixPos < 10) {
      features = pPSM->features;
      for (unsigned int a = 0; a < numFeatures; ++a) {
        cerr << features[a] << " ";
      }
      cerr << endl;
    }

  }
  for (it = setVec.begin(); it != setVec.end(); ++it) {
    int ixPos = -1;
    while ((pPSM = (*it)->getNext(ixPos)) != NULL) {
      features = pPSM->features;
      normalize(features, features, 0, numFeatures);
      features = pPSM->retentionFeatures;
      normalize(features, features, numFeatures, numRetentionFeatures);
    }
  }
  if (VERB > 4) {
    cerr << "First 10 feature vectors after normalization" << endl;
    cerr.precision(3);
    it = setVec.begin();
    int ixPos = -1;
    while ((pPSM = (*it)->getNext(ixPos)) != NULL && ixPos < 10) {
      features = pPSM->features;
      for (unsigned int a = 0; a < numFeatures; ++a) {
        cerr << features[a] << " ";
      }
      cerr << endl;
    }

  }
}

void Normalizer::normalize(const double *in, double* out, size_t offset,
                           size_t numFeatures) {
  for (unsigned int ix = 0; ix < numFeatures; ++ix) {
    out[ix] = (in[ix] - sub[offset + ix]) / div[offset + ix];
  }
}

int Normalizer::subclass_type = STDV;
Normalizer* Normalizer::theNormalizer = NULL;

Normalizer * Normalizer::getNormalizer() {
  if (theNormalizer == NULL) {
    if (subclass_type == UNI)
      theNormalizer = new UniNormalizer();
    else
      theNormalizer = new StdvNormalizer();
  }
  return theNormalizer;
}

void Normalizer::setType(int type) {
  assert(type == UNI || type == STDV);
  subclass_type = type;
}
