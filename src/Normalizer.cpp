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
#include <assert.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
using namespace std;
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "UniNormalizer.h"
#include "Globals.h"

int Normalizer::subclass_type = STDV;
Normalizer* Normalizer::theNormalizer = NULL;

Normalizer::Normalizer() {
}

Normalizer::~Normalizer() {
}

void Normalizer::normalizeSet(vector<double*>& featuresV,
                              vector<double*>& rtFeaturesV) {
  normalizeSet(featuresV, 0, numFeatures);
  normalizeSet(rtFeaturesV, numFeatures, numRetentionFeatures);
}

void Normalizer::normalizeSet(vector<double*>& featuresV,
                              size_t offset, size_t numFeatures) {
  double* features;
  vector<double*>::iterator it = featuresV.begin();
  for (; it != featuresV.end(); ++it) {
    features = *it;
    normalize(features, features, offset, numFeatures);
  }
}

void Normalizer::normalize(const double* in, double* out, size_t offset,
                           size_t numFeatures) {
  for (unsigned int ix = 0; ix < numFeatures; ++ix) {
    out[ix] = (in[ix] - sub[offset + ix]) / div[offset + ix];
  }
}

Normalizer* Normalizer::getNormalizer() {
  if (theNormalizer == NULL) {
    if (subclass_type == UNI) {
      theNormalizer = new UniNormalizer();
    } else {
      theNormalizer = new StdvNormalizer();
    }
  }
  return theNormalizer;
}

void Normalizer::setType(int type) {
  assert(type == UNI || type == STDV);
  subclass_type = type;
}
