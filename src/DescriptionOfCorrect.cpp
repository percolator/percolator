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
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include "Globals.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "DescriptionOfCorrect.h"
#include "MassHandler.h"
#include "PSMDescription.h"
#include "Enzyme.h"

string DescriptionOfCorrect::isoAlphabet = "DECYHKR";
unsigned int DescriptionOfCorrect::docFeatures = 15;
float DescriptionOfCorrect::pKiso[7] = { -3.86f, -4.25f, -8.33f, -10.0f, 6.0f,
                                         10.5f, 12.4f }; // Lehninger
float DescriptionOfCorrect::pKN = 9.69f;
float DescriptionOfCorrect::pKC = 2.34f;

DescriptionOfCorrect::DescriptionOfCorrect() {
}

DescriptionOfCorrect::~DescriptionOfCorrect() {
}

void DescriptionOfCorrect::calcRegressionFeature(PSMDescription* psm) {
  string peptide = psm->getFullPeptide();
  string::size_type pos1 = peptide.find('.');
  string::size_type pos2 = peptide.find('.', ++pos1);
  string pep = peptide.substr(pos1, pos2 - pos1);
  psm->setIsoElectricPoint(isoElectricPoint(pep));
  if (psm->getRetentionFeatures()) {
    RTModel::fillFeaturesAllIndex(pep, psm->getRetentionFeatures());
  }
  //cout <<  peptide << " " << pep << " " << psm->getRetentionFeatures()[0] << endl;
}

void DescriptionOfCorrect::trainCorrect() {
  // Get rid of redundant peptides
  sort(psms.begin(), psms.end(), PSMDescription::ptrLess);
  psms.erase(std::unique(psms.begin(), psms.end(), PSMDescription::ptrEqual), psms.end());
  
  // Get averages
  double piSum = 0.0, dMSum = 0.0;
  for (size_t ix = 0; ix < psms.size(); ++ix) {
    piSum += psms[ix]->getIsoElectricPoint();
    dMSum += psms[ix]->getMassDiff();
  }
  if (psms.size() == 0) {
    avgPI = 0.0;
    avgDM = 0.0;
  } else {
    avgPI = piSum / psms.size();
    avgDM = dMSum / psms.size();
  }
  rtModel.trainRetention(psms);
  if (VERB > 2) {
    cerr << "Description of correct recalibrated, avg pI=" << avgPI
        << " avg dM=" << avgDM << endl;
  }
}

void DescriptionOfCorrect::setFeatures(PSMDescription* psm) {
  assert(DataSet::getFeatureNames().getDocFeatNum() > 0);
  psm->setPredictedRetentionTime(rtModel.estimateRT(psm->getRetentionFeatures()));
  size_t docFeatNum = DataSet::getFeatureNames().getDocFeatNum();
  double dm = abs(psm->getMassDiff() - avgDM);
  double drt = abs(psm->getRetentionTime() - psm->getPredictedRetentionTime());
  if (docFeatures & 1) {
    psm->features[docFeatNum] = abs(psm->getIsoElectricPoint() - avgPI);
  }
  if (docFeatures & 2) {
    psm->features[docFeatNum + 1] = dm;
  }
  if (docFeatures & 4) {
    psm->features[docFeatNum + 2] = drt;
  }
  // double ddrt=drt/(1+log(max(1.0,PSMDescriptionDOC::unnormalize(psm->getRetentionTime))));
  if (docFeatures & 8) {
    psm->features[docFeatNum + 3] = sqrt(dm * drt);
  }
}

void DescriptionOfCorrect::setFeaturesNormalized(PSMDescription* psm, Normalizer* pNorm) {
  setFeatures(psm);
  size_t docFeatNum = DataSet::getFeatureNames().getDocFeatNum();
  if (docFeatures & 1) {
    psm->features[docFeatNum] = pNorm->normalize(psm->features[docFeatNum], docFeatNum);
  }
  if (docFeatures & 2) {
    psm->features[docFeatNum + 1] = pNorm->normalize(psm->features[docFeatNum + 1], docFeatNum + 1);
  }
  if (docFeatures & 4) {
    psm->features[docFeatNum + 2] = pNorm->normalize(psm->features[docFeatNum + 2], docFeatNum + 2);
  }
  // double ddrt=drt/(1+log(max(1.0,PSMDescription::unnormalize(psm->getRetentionTime))));
  if (docFeatures & 8) {
    psm->features[docFeatNum + 3] = pNorm->normalize(psm->features[docFeatNum + 3], docFeatNum + 3);
  }
}

double DescriptionOfCorrect::isoElectricPoint(const string& pep) {
  // Overall amino acid composition features
  string::size_type pos = isoAlphabet.size();
  vector<unsigned int> numAA(pos);
  for (string::const_iterator it = pep.begin(); it != pep.end(); it++) {
    pos = isoAlphabet.find(*it);
    if (pos != string::npos) {
      numAA[pos]++;
    }
  }
  double pH = 6.5, pHlow = 2.0, pHhigh = 13.0;
  double epsilon = 0.01;
  while ((pH - pHlow > epsilon) || (pHhigh - pH > epsilon)) {
    double NQ = 1 / (1 + pow(10, (pH - pKN))) - 1 / (1 + pow(10,
                                                             (pKC - pH)));
    for (size_t ix = 0; ix < numAA.size(); ix++) {
      if (numAA[ix] == 0) {
        continue;
      }
      if (pKiso[ix] > 0) {
        NQ += numAA[ix] / (1 + pow(10, (pH - pKiso[ix])));
      } else {
        NQ -= numAA[ix] / (1 + pow(10, (-pKiso[ix] - pH)));
      }
    }
    if (NQ < 0) { //Bisection method
      pHhigh = pH;
      pH -= ((pH - pHlow) / 2);
    } else {
      pHlow = pH;
      pH += ((pHhigh - pH) / 2);
    }
  }
  return pH;
}
