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
#include <cstring>
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
float DescriptionOfCorrect::pKiso[7] = {-3.86,-4.25,-8.33,-10.0,6.0,10.5,12.4}; // Lehninger
float DescriptionOfCorrect::pKN = 9.69;
float DescriptionOfCorrect::pKC = 2.34;

DescriptionOfCorrect::DescriptionOfCorrect()
{
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}

void DescriptionOfCorrect::calcRegressionFeature(PSMDescription &psm) {
  string peptide = psm.getFullPeptide();
  string::size_type pos1 = peptide.find('.');
  string::size_type pos2 = peptide.find('.',++pos1);
  string pep = peptide.substr(pos1,pos2-pos1);
  psm.pI = isoElectricPoint(pep);
  if (psm.getRetentionFeatures()) {
    RTModel::fillFeaturesAllIndex(pep, psm.getRetentionFeatures());
  }
//  cout <<  peptide << " " << pep << " " << retentionFeatures[0] << endl;
}

void DescriptionOfCorrect::trainCorrect() {
  // Get rid of redundant peptides
  sort(psms.begin(),psms.end(),less<PSMDescription>());
  psms.resize(distance(psms.begin(),unique(psms.begin(),psms.end())));

  // Get rid of non enzymatic peptides

  remove_if(psms.begin(),psms.end(),mem_fun_ref(&PSMDescription::isNotEnzymatic));

  // Get averages
  double piSum=0.0, dMSum=0.0;
  for(size_t ix=0; ix<psms.size(); ++ix) {
    piSum += psms[ix].pI;
    dMSum += psms[ix].massDiff;
  }
  if (psms.size()==0) {
    avgPI = 0.0; avgDM = 0.0;
  } else {
    avgPI = piSum/psms.size();
    avgDM = dMSum/psms.size();
  }
  rtModel.trainRetention(psms);

  if (VERB>2) cerr << "Description of correct recalibrated, avg pI=" << avgPI << " avg dM=" << avgDM << endl;

}
void DescriptionOfCorrect::setFeatures(PSMDescription& psm) {
  assert(DataSet::getFeatureNames().getDocFeatNum()>0);
  psm.predictedTime=rtModel.estimateRT(psm.retentionFeatures);
  size_t docFeatNum = DataSet::getFeatureNames().getDocFeatNum();
  double dm = abs(psm.massDiff-avgDM);
  double drt=abs(psm.retentionTime-psm.predictedTime);

  if (docFeatures & 1) psm.features[docFeatNum] = Normalizer::getNormalizer()->normalize(abs(psm.pI-avgPI),docFeatNum);
  else psm.features[docFeatNum] = 0.0;

  if (docFeatures & 2) psm.features[docFeatNum+1] = Normalizer::getNormalizer()->normalize(dm,docFeatNum+1);
  else psm.features[docFeatNum+1] = 0.0;


  if (docFeatures & 4) psm.features[docFeatNum+2] = Normalizer::getNormalizer()->normalize(drt,docFeatNum+2);
  else psm.features[docFeatNum+2] = 0.0;

  // double ddrt=drt/(1+log(max(1.0,PSMDescription::unnormalize(psm.retentionTime))));
  if (docFeatures & 8) psm.features[docFeatNum+3] = Normalizer::getNormalizer()->normalize(sqrt(dm*drt),docFeatNum+3);
  else psm.features[docFeatNum+3] = 0.0;

}

double DescriptionOfCorrect::isoElectricPoint(const string& pep) {
  // Overall amino acid composition features
  string::size_type pos = isoAlphabet.size();
  vector<unsigned int> numAA(pos);
  for (string::const_iterator it=pep.begin();it!=pep.end();it++) {
    pos=isoAlphabet.find(*it);
    if (pos!=string::npos) numAA[pos]++;
  }
  double pH = 6.5, pHlow = 2.0, pHhigh = 13.0;
  double epsilon = 0.01;

  while((pH-pHlow > epsilon) || (pHhigh-pH > epsilon)) {
    double NQ = 1/(1+pow(10,(pH-pKN))) - 1/(1+pow(10,(pKC-pH)));
    for(size_t ix=0; ix<numAA.size();ix++) {
      if (numAA[ix]==0)
        continue;
      if (pKiso[ix]>0) {
        NQ += numAA[ix]/(1+pow(10,(pH-pKiso[ix])));
      } else {
        NQ -= numAA[ix]/(1+pow(10,(-pKiso[ix]-pH)));
      }
    }
    if(NQ<0) {  //Bisection method
        pHhigh = pH;
        pH -= ((pH-pHlow)/2);
    } else {
        pHlow = pH;
        pH += ((pHhigh-pH)/2);
    }
 }
 return pH;
}
/*
<<<<<<< HEAD:src/DescriptionOfCorrect.cpp
=======

void DescriptionOfCorrect::fillFeaturesAllIndex(const string& pep, double *features) {
//  unsigned int ptms=0;
  string peptide = pep;
  string::size_type posP = 0;
  vector<unsigned int> ptms(ptmAlphabet.length(),0u);
  while (posP < ptmAlphabet.length()) {
    string::size_type pos = 0;
    while ( (pos = peptide.find(ptmAlphabet[posP], pos)) != string::npos ) {
      peptide.replace( pos, 1, "");
      ++pos;
      ++ptms[posP];
    }
    ++posP;
  }

  if(!doKlammer) {
    features = fillFeaturesIndex(peptide, krokhin_index, features);
    features = fillFeaturesIndex(peptide, krokhin100_index, features);
    features = fillFeaturesIndex(peptide, krokhinTFA_index, features);
//    features = fillFeaturesIndex(peptide, hessa_index, features);
//    features = fillFeaturesIndex(peptide, kytedoolittle_index, features);
    *(features++) = peptide.size();
    for(vector<unsigned int>::iterator ptm=ptms.begin();ptm != ptms.end();++ptm)
       *(features++) = (double) *ptm;
    features = fillAAFeatures(peptide, features);
  } else {
  	// Klammer et al. features
    features = fillAAFeatures(peptide, features);
    features = fillAAFeatures(peptide.substr(0,1), features);
    features = fillAAFeatures(peptide.substr(peptide.size()-2,1), features);
    char Ct = peptide[peptide.size()-1];
    *(features++) = Enzyme::isEnzymatic(Ct,'A');
    *(features++) = peptide.size();
    *(features++) = indexSum(aa_weights,peptide)+1.0079+17.0073; //MV

  }
//  cout <<  pep << " " << peptide << endl;
}


float DescriptionOfCorrect::krokhin_index['Z'-'A'+1] =
         {1.1, 0.0, 0.45, 0.15, 0.95,  10.9, -0.35, -1.45, 8.0, 0.0, -2.05, 9.3, 6.2,
          -0.85,0.0,2.1,-0.4,-1.4,-0.15,0.65,0.0,5.0,12.25,0.0,4.85,0.0};
float DescriptionOfCorrect::krokhin100_index['Z'-'A'+1] =
         {1,0,0.1,0.15,1,11.67,-0.35,-3.0,7.95,0,-3.4,9.4,6.25,
         -0.95,0,1.85,-0.6,-2.55,-0.15,0.65,0,4.7,13.35,0,5.35,0};
float DescriptionOfCorrect::krokhinC2_index['Z'-'A'+1] =
         {0.5,0.0,0.2,0.4,0,9.5,0.15,-0.2,6.6,0.0,-1.5,7.4,5.7,-0.2,0.0,
          2.1,-0.2,-1.1,-0.1,0.6,0.0,3.4,11.8,0.0,4.5,0.0};
float DescriptionOfCorrect::krokhinTFA_index['Z'-'A'+1] =
         {1.11,0.0,0.04,-0.22,1.08,11.34,-0.35,-3.04,7.86,0.0,-3.53,9.44,6.57,
          -1.44,0.0,1.62,-0.53,-2.58,-0.33,0.48,0.0,4.86,13.12,0.0,5.4,0.0};

// negated hessa scale
float DescriptionOfCorrect::hessa_index['Z'-'A'+1] =
         {-0.11,-0.0,0.13,-3.49,-2.68,0.32,-0.74,-2.06,0.60,0.0,-2.71,0.55,0.10,-2.05,
          0.0,-2.23,-2.36,-2.58,-0.84,-0.52,0.0,0.31,-0.30,0.0,-0.68,0.0};
float DescriptionOfCorrect::kytedoolittle_index['Z'-'A'+1] =
         {1.80,0.0,2.50,-3.50,-3.50,2.80,-0.40,-3.20,4.50,0.0,-3.90,3.80,1.90,-3.50,
          0.0,-1.60,-3.50,-4.50,-0.80,-0.70,0.0,4.20,-0.90,0.0,-1.30,0.0};
float DescriptionOfCorrect::aa_weights['Z'-'A'+1] = {
  71.0788,0, 103.1448,115.0886,129.1155,147.1766,57.052,137.1412,113.1595,0,128.1742,113.1595,131.1986,114.1039,0,97.1167,128.1308,156.1876,87.0782,101.1051,0,99.1326,186.2133,0,163.1760,0};
>>>>>>> bf133cae15dafdae12f5da41f6a0f92a2d28e852:src/DescriptionOfCorrect.cpp
*/
