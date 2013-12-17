
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

#include "FeatureNames.h"
#include "Globals.h"

size_t FeatureNames::numFeatures = 0;

FeatureNames::FeatureNames() {
  minCharge = 100;
  maxCharge = -1;
  chargeFeatNum = -1;
  enzFeatNum = -1;
  numSPFeatNum = -1;
  ptmFeatNum = -1;
  intraSetFeatNum = -1;
  quadraticFeatNum = -1;
  docFeatNum = -1;
}

FeatureNames::~FeatureNames() {
}

void FeatureNames::setFromXml( const ::percolatorInNs::featureDescriptions & fdes, bool calcDOC ) {
  //assert(featureNames.empty());
  for( const auto & descr : fdes.featureDescription() ) {
    featureNames.push_back(descr.name());
  }

  if (calcDOC) {
    docFeatNum = featureNames.size();
    featureNames.push_back("docpI");
    featureNames.push_back("docdM");
    featureNames.push_back("docRT");
    featureNames.push_back("docdMdRT");
  }
  setNumFeatures(featureNames.size());
  if (VERB>2) {
    std::cerr << "in FeatureNames::setFromXml\n";
  }
  if (VERB>1) {
    std::cerr << "Features:\n";
    std::copy( featureNames.begin(), featureNames.end(),
        std::ostream_iterator<std::string>(std::cerr, " "));
    std::cerr << "\n";
  }
  if (VERB>2) {
    std::cerr << "end of FeatureNames::setFromXml\n";
  }
  return;
}

string FeatureNames::getFeatureNames(bool skipDOC) {
  int n = (skipDOC && docFeatNum > 0) ? docFeatNum
      : (int)featureNames.size();
  ostringstream oss;
  if (!featureNames.empty()) {
    int featNum = 0;
    oss << featureNames[featNum++];
    for (; featNum < n; ++featNum) {
      oss << "\t" << featureNames[featNum];
    }
  }
  return oss.str();
}

void FeatureNames::setFeatures(string& line, size_t skip, size_t numFields) {
  if (!featureNames.empty()) {
    return;
  }
  istringstream iss(line);
  string tmp;
  while (iss.good() && skip && --skip) {
    iss >> tmp;
  }
  int remain = numFields;
  while (iss.good() && remain && --remain) {
    iss >> tmp;
    featureNames.push_back(tmp);
  }
  assert(featureNames.size() == numFields);
  setNumFeatures(featureNames.size());
}
