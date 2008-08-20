#include <sstream>
#include "FeatureNames.h"


size_t FeatureNames::numFeatures = 0;


FeatureNames::FeatureNames()
{
  minCharge = 100;
  maxCharge = -1;
  chargeFeatNum = -1;
  enzFeatNum = -1;
  numSPFeatNum = -1;
  ptmFeatNum = -1;
  rank1FeatNum = -1;
  aaFeatNum = -1;
  intraSetFeatNum = -1;
  quadraticFeatNum = -1;
  docFeatNum = -1;

}

FeatureNames::~FeatureNames()
{
}

string FeatureNames::getFeatureNames() {
  ostringstream oss;
  if (!featureNames.empty()) {
    vector<string>::iterator featNum = featureNames.begin(); 
    oss << *(featNum++);
    for (; featNum!=featureNames.end(); ++featNum )
      oss << "\t" << *featNum;
  }
  return oss.str();
}


void FeatureNames::setSQTFeatures(
  int minC, int maxC, 
  bool doEnzyme, 
  bool calcPTMs, 
  bool doManyHitsPerSpectrum, 
  const string& aaAlphabet, 
  bool calcQuadratic, 
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
  minCharge = minC; maxCharge = maxC; 
  for(int charge=minCharge; charge <= maxCharge; ++charge) {
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
  if (doManyHitsPerSpectrum) {
    rank1FeatNum = featureNames.size();
    featureNames.push_back("rank1");
  }
  if (!aaAlphabet.empty()) {
    aaFeatNum = featureNames.size();
    for (string::const_iterator it=aaAlphabet.begin();it!=aaAlphabet.end();it++)
      featureNames.push_back(*it + "-Freq");
  }
  if(calcQuadratic) {
    quadraticFeatNum = featureNames.size();
    for(int f1=1;f1<quadraticFeatNum;++f1) {
      for(int f2=0;f2<f1;++f2) {
        ostringstream feat;
        feat << "f" << f1+1 << "*" << "f" << f2+1;
        featureNames.push_back(feat.str());
      }    
    }
  }
  if (calcDOC) {
    docFeatNum = featureNames.size();
    featureNames.push_back("docpI");
    featureNames.push_back("docdM");
    featureNames.push_back("docRT");
  }
  setNumFeatures(featureNames.size());  
}
