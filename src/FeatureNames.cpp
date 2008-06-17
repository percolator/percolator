#include <sstream>
#include "FeatureNames.h"

FeatureNames::FeatureNames()
{
  minCharge = -1;
  maxCharge = -1;
  chargeFeatNum = -1;
  enzFeatNum = -1;
  numSPFeatNum = -1;
  ptmFeatNum = -1;
  rank1FeatNum = -1;
  aaFeatNum = -1;
  intraSetFeatNum = -1;
  docFeatNum = -1;

}

FeatureNames::~FeatureNames()
{
}

string getFeatureNames() {
  ostringstream oss;
  if (!featureNames.empty()) {
    vector<string>::iterator featNam = featureNames.begin(); 
    oss << *(featNam++);
    for (; featNum!=featureNames.end(); ++featNum )
      oss << "\t" << *(featNam++);
  }
  return oss.str();
}


string DataSet::setSQTFeatures(
  int minCharge, int maxCharge, 
  bool doEnzyme, 
  bool calcPTMs, 
  bool doManyHitsPerSpectrum, 
  bool calcAAFrequencies, 
  bool calcIntraSetFeatures, 
  bool calcDOC)
{
  if (!featureNames.empty())
    return;
  featureNames.push_back("lnrSp");
  featureNames.push_back("deltLCn");
  featureNames.push_back("deltCn");
  featureNames.push_back("Xcorr");
  featureNames.push_back("Sp");
  featureNames.push_back("IonFrac");
  featureNames.push_back("Mass");
  featureNames.push_back("PepLen");
  chargeFeatNum = featureNames.size();
  for(int charge=minCharge; charge <= maxCharge; ++charge)
    featureNames.push_back("Charge" + charge);
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
  if (doManyHitsPerSpectrum) {
    rank1FeatNum = featureNames.size();
    featureNames.push_back("rank1");
  }
  if (calcAAFrequencies) {
    aaFeatNum = featureNames.size();
    for (string::const_iterator it=aaAlphabet.begin();it!=aaAlphabet.end();it++)
      featureNames.push_back(*it + "-Freq");
  }
  if (calcIntraSetFeatures) {
    intraSetFeatNum = featureNames.size();
    featureNames.push_back("numPep");
    featureNames.push_back("pepSite");
  }
  if (calcDOC) {
    docFeatNum = featureNames.size();
    featureNames.push_back("docpI");
    featureNames.push_back("docdM");
    featureNames.push_back("docRT");
  }  
}
