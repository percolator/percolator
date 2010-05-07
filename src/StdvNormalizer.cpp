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
#include <vector>
#include <iostream>
#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif
#include <math.h>
#include <set>
#include <vector>
#include <string>
using namespace std;
#include "PSMDescription.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "SetHandler.h"
#include "Globals.h"

StdvNormalizer::StdvNormalizer() {
}

StdvNormalizer::~StdvNormalizer() {
}

void StdvNormalizer::unnormalizeweight(const vector<double>& in, vector<
    double>& out) {
  double sum = 0;
  unsigned int i = 0;
  for (; i < FeatureNames::getNumFeatures(); i++) {
    out[i] = in[i] / div[i];
    sum += sub[i] * in[i] / div[i];
  }
  out[i] = in[i] - sum;
}

void StdvNormalizer::normalizeweight(const vector<double>& in, vector<
    double>& out) {
  double sum = 0;
  size_t i = 0;
  for (; i < FeatureNames::getNumFeatures(); i++) {
    out[i] = in[i] * div[i];
    sum += sub[i] * in[i];
  }
  out[i] = in[i] + sum;
}

void StdvNormalizer::setSet(set<DataSet *> & setVec, size_t nf, size_t nrf) {
  numFeatures = nf;
  numRetentionFeatures = nrf;
  sub.resize(nf + nrf, 0.0);
  div.resize(nf + nrf, 0.0);
  double n = 0.0;
  double * features;
  PSMDescription* pPSM;
  size_t ix;
  set<DataSet *>::iterator it;
  for (it = setVec.begin(); it != setVec.end(); ++it) {
    assert(*it);
    int ixPos = -1;
    while ((pPSM = (*it)->getNext(ixPos)) != NULL) {
      assert(pPSM);
      n++;
      features = pPSM->features;
      for (ix = 0; ix < numFeatures; ++ix) {
        sub[ix] += features[ix];
      }
      features = pPSM->retentionFeatures;
      for (; ix < numFeatures + numRetentionFeatures; ++ix) {
        sub[ix] += features[ix - numFeatures];
      }
    }
  }
  if (VERB > 2) {
    cerr.precision(2);
    cerr << "Normalization factors" << endl << "Type\t"
        << DataSet::getFeatureNames().getFeatureNames() << endl << "Avg ";
  }
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (n > 0.0) sub[ix] /= n;
    if (VERB > 2) cerr << "\t" << sub[ix];
  }
  for (it = setVec.begin(); it != setVec.end(); ++it) {
    int ixPos = -1;
    while ((pPSM = (*it)->getNext(ixPos)) != NULL) {
      features = pPSM->features;
      for (ix = 0; ix < numFeatures; ++ix) {
        if (!isfinite(features[ix])) cerr
            << "Reached strange feature with val=" << features[ix]
            << " at row=" << ix << ", col=" << ixPos << endl;
        double d = features[ix] - sub[ix];
        div[ix] += d * d;
      }
      features = pPSM->retentionFeatures;
      for (; ix < numFeatures + numRetentionFeatures; ++ix) {
        if (!isfinite(features[ix - numFeatures])) cerr
            << "Reached strange feature with val=" << features[ix
            - numFeatures] << " at row=" << ix << ", col=" << ixPos
            << endl;
        double d = features[ix - numFeatures] - sub[ix];
        div[ix] += d * d;
      }
    }
  }
  if (VERB > 2) cerr << endl << "Stdv";
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (div[ix] <= 0 || n == 0) {
      div[ix] = 1.0;
    } else {
      div[ix] = sqrt(div[ix] / n);
    }
    if (VERB > 2) cerr << "\t" << div[ix];
  }
  if (VERB > 2) cerr << endl;
}
