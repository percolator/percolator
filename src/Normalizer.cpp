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

void Normalizer::normalizeSet(vector<double*> & featuresV,
                              vector<double*> & rtFeaturesV) {
  double* features;
  vector<double*>::iterator it = featuresV.begin();
  for (; it != featuresV.end(); ++it) {
    features = *it;
    normalize(features, features, 0, numFeatures);
  }
  vector<double*>::iterator rtit = rtFeaturesV.begin();
  for (; rtit != rtFeaturesV.end(); ++rtit) {
    features = *rtit;
    normalize(features, features, numFeatures, numRetentionFeatures);
  }
}

void Normalizer::normalize(const double* in, double* out, size_t offset,
                           size_t numFeatures) {
  for (unsigned int ix = 0; ix < numFeatures; ++ix) {
    out[ix] = (in[ix] - sub[offset + ix]) / div[offset + ix];
  }
}

void Normalizer::unNormalizeSet(vector<double*> & rtFeaturesV) {
  double* features;
  for (int i = 0; i < rtFeaturesV.size(); ++i) {
    features = rtFeaturesV[i];
    for (int j = 0; j < numRetentionFeatures; ++j) {
      features[j] = (features[j] * div[j]) + sub[j];
    }
  }
}
// normalize a set of PSMs
/*
 void Normalizer::normalizeSet(vector<PSMDescription> & psms)
 {
 vector<PSMDescription>::iterator it;
 double * retFeatures;

 cout << "Normalizing..." << endl;
 for(it = psms.begin(); it != psms.end(); ++it)
 {
 retFeatures = it->retentionFeatures;
 normalize(retFeatures, retFeatures, 0, numRetentionFeatures);
 }
 cout << "Done." << endl << endl;
 }
 */
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
