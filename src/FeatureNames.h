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
#ifndef FEATURENAMES_H_
#define FEATURENAMES_H_

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <boost/foreach.hpp>

#include "percolator_in.hxx"

using namespace std;

class FeatureNames {
  public:
    FeatureNames();
    virtual ~FeatureNames();
    string getFeatureNames(bool skipDOC = false);
    static inline size_t getNumFeatures() {
      return numFeatures;
    }
    static inline void setNumFeatures(size_t nf) {
      if (!numFeatures) {
        numFeatures = nf;
      }
    }


    void setFromXml( const ::percolatorInNs::featureDescriptions & fdes, bool calcDOC );

    void insertFeature(const string& featureName) {
      if(std::find(featureNames.begin(), featureNames.end(), featureName)==featureNames.end())
	featureNames.push_back(featureName);
    }
    void setFeatures(string& line, size_t skip, size_t numFeatures);

    int getFeatureNumber(const string& featureName) {
      return distance(find(featureNames.begin(),
                           featureNames.end(),
                           featureName), featureNames.begin());
    }
    // SQT Feature Number getters, will return -1 if not defined.
    int getMinCharge() {
      return minCharge;
    }
    int getMaxCharge() {
      return maxCharge;
    }
    int getChargeFeatNum() {
      assert(1==0);
      return chargeFeatNum;
    }
    int getEnzFeatNum() {
      assert(1==0);
      return enzFeatNum;
    }
    int getNumSPFeatNum() {
      assert(1==0);
      return numSPFeatNum;
    }
    int getPtmFeatNum() {
      assert(1==0);
      return ptmFeatNum;
    }
    int getIntraSetFeatNum() {
      assert(1==0);
      return intraSetFeatNum;
    }
    int getQuadraticFeatNum() {
      assert(1==0);
      return quadraticFeatNum;
    }
    int getDocFeatNum() {
      return docFeatNum;
    }
    void setDocFeatNum(int fn) {
      docFeatNum = fn;
    }
  protected:
    vector<string> featureNames;
    static size_t numFeatures;
    int minCharge, maxCharge;
    int chargeFeatNum, enzFeatNum, numSPFeatNum, ptmFeatNum,
        intraSetFeatNum, quadraticFeatNum, docFeatNum;
};

#endif /*FEATURENAMES_H_*/
