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
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"
#include "SqtSanityCheck.h"

// Class doing the sanity check for non SQT file condition
// In most sence a place holder as very little logic is build-in in this case

SanityCheck::SanityCheck() :
  initPositives(0), pTestset(NULL), pTrainset(NULL) {
}

SanityCheck::~SanityCheck() {
}

bool SanityCheck::overRule = false;
string SanityCheck::initWeightFN = "";
int SanityCheck::initDefaultDir = 0;

/**
 * Returns an instance of the appropriate sanity check based on information
 * contained in otherCall.
 */
SanityCheck* SanityCheck::initialize(string otherCall){
  if(initWeightFN!="" || initDefaultDir!=0) {
    return new SanityCheck();
  }
  else if (otherCall.find(SqtSanityCheck::fingerPrint)!= string::npos){
    return new SqtSanityCheck();
  } else return new SanityCheck();
}

int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr) {
  pTestset = &testset;
  pTrainset = &trainset;
  fdr = test_fdr;
  if (initWeightFN.size() > 0) {
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN.data(), ios::in);
    readWeights(weightStream, ww);
    weightStream.close();
    assert(pNorm);
    pNorm->normalizeweight(ww, w[0]);
    for (size_t set = 1; set < w.size(); ++set) {
      copy(w[0].begin(), w[0].end(), w[set].begin());
    }
  } else {
    getDefaultDirection(w);
  }
  initPositives = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    initPositives
        += (*pTrainset)[set].getInitDirection(fdr, w[set], false);
  }
  initPositives /= (int)(max<double> ((double)(w.size() - 1), 1.0));
  return initPositives;
}

void SanityCheck::getDefaultDirection(vector<vector<double> >& w) {
  if (!initDefaultDir) {
    // Set init direction to be the most discriminative direction
    for (size_t set = 0; set < w.size(); ++set) {
      (*pTrainset)[set].getInitDirection(fdr, w[set], true);
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      for (size_t ix = 0; ix < w[set].size(); ix++) {
        w[set][ix] = 0;
      }
      w[set][abs(initDefaultDir) - 1] = (initDefaultDir < 0 ? -1 : 1);
    }
  }
}

bool SanityCheck::validateDirection(vector<vector<double> >& w) {
  if (!pTestset) {
    cerr << "SanityCheck wrongly configured" << endl;
    exit(-1);
  }
  int overFDR = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    overFDR += (*pTestset)[set].calcScores(w[set], fdr);
  }
  if (VERB > 0) {
    cerr << "Found " << overFDR << " target PSMs scoring over " << fdr
        * 100 << "% FDR level on testset" << endl;
  }
  if (overFDR <= 0) {
    cerr << "No target score better than best decoy" << endl;
    resetDirection(w);
    return false;
  }
  if (initPositives >= overFDR) {
    cerr << "Less identifications (" << overFDR << " vs " << initPositives
        << ") after percolator processing than before processing" << endl;
    resetDirection(w);
    return false;
  }
  if (initDefaultDir) {
    for (size_t set = 0; set < w.size(); ++set) {
      if (w[set][abs(initDefaultDir) - 1] * initDefaultDir <= 0) {
        if (VERB > 1) {
          cerr
            << "Warning, wrong sign of the weight for main scoring direction"
            << endl;
        }
        resetDirection(w);
        return false;
      }
    }
  }
  return true;
}

void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  char buffer[1024], c;
  while (!(((c = weightStream.get()) == '-') || (c >= '0' && c <= '9'))) {
    weightStream.getline(buffer, 1024);
  }
  weightStream.getline(buffer, 1024);
  // Get second line containing raw features
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream >> w[ix];
  }
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<vector<double> >& w) {
  if (!overRule) {
    cerr << "Reseting score vector, using default vector" << endl;
    getDefaultDirection(w);
  }
}

