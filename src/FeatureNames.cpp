
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

#include "FeatureNames.h"

size_t FeatureNames::numFeatures = 0;

FeatureNames::FeatureNames() {
  minCharge = 100;
  maxCharge = -1;
  chargeFeatNum = -1;
  enzFeatNum = -1;
  numSPFeatNum = -1;
  ptmFeatNum = -1;
  pngFeatNum = -1;
  aaFeatNum = -1;
  intraSetFeatNum = -1;
  quadraticFeatNum = -1;
  docFeatNum = -1;
}

FeatureNames::~FeatureNames() {
}

void FeatureNames::setFromXml( const ::percolatorInNs::featureDescriptions & fdes, bool calcDOC ) {
  //assert(featureNames.empty());
  BOOST_FOREACH( const ::percolatorInNs::featureDescription & descr,  fdes.featureDescription() ) {
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
  std::cout << "in FeatureNames::setFromXml" << std::endl;
  std::copy( featureNames.begin(), featureNames.end(), std::ostream_iterator<std::string>(std::cout, " "));
  std::cout << "end of FeatureNames::setFromXml" << std::endl;
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

void FeatureNames::setSQTFeatures(int minC, int maxC, bool doEnzyme,
                                  bool calcPTMs, bool doPNGaseF,
                                  const string& aaAlphabet,
                                  bool calcQuadratic, bool calcDOC) {
  if (!featureNames.empty()) {
    return;
  }
  featureNames.push_back("lnrSp");
  featureNames.push_back("deltLCn");
  featureNames.push_back("deltCn");
  featureNames.push_back("Xcorr");
  featureNames.push_back("Sp");
  featureNames.push_back("IonFrac");
  featureNames.push_back("Mass");
  featureNames.push_back("PepLen");
  chargeFeatNum = featureNames.size();
  minCharge = minC;
  maxCharge = maxC;
  for (int charge = minCharge; charge <= maxCharge; ++charge) {
    ostringstream cname;
    cname << "Charge" << charge;
    featureNames.push_back(cname.str());
  }
  if (doEnzyme) {
    enzFeatNum = featureNames.size();
    featureNames.push_back("enzN");
    featureNames.push_back("enzC");
    featureNames.push_back("enzInt");
  }
  numSPFeatNum = featureNames.size();
  featureNames.push_back("lnNumSP");
  featureNames.push_back("dM");
  featureNames.push_back("absdM");
  if (calcPTMs) {
    ptmFeatNum = featureNames.size();
    featureNames.push_back("ptm");
  }
  if (doPNGaseF) {
    pngFeatNum = featureNames.size();
    featureNames.push_back("PNGaseF");
  }
  if (!aaAlphabet.empty()) {
    aaFeatNum = featureNames.size();
    for (string::const_iterator it = aaAlphabet.begin(); it
        != aaAlphabet.end(); it++) {
      featureNames.push_back(*it + "-Freq");
    }
  }
  if (calcQuadratic) {
    quadraticFeatNum = featureNames.size();
    for (int f1 = 1; f1 < quadraticFeatNum; ++f1) {
      for (int f2 = 0; f2 < f1; ++f2) {
        ostringstream feat;
        feat << "f" << f1 + 1 << "*" << "f" << f2 + 1;
        featureNames.push_back(feat.str());
      }
    }
  }
  if (calcDOC) {
    docFeatNum = featureNames.size();
    featureNames.push_back("docpI");
    featureNames.push_back("docdM");
    featureNames.push_back("docRT");
    featureNames.push_back("docdMdRT");
  }
  setNumFeatures(featureNames.size());
  std::cout << "in FeatureNames::setSQTFeatures" << std::endl;
  std::copy( featureNames.begin(), featureNames.end(), std::ostream_iterator<std::string>(std::cout, " "));
  std::cout << "end of FeatureNames::setSQTFeatures" << std::endl;
  return;
}
