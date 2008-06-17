#include<vector>
#include<cmath>
#include "DataSet.h"
#include "DescriptionOfCorrect.h"



DescriptionOfCorrect::DescriptionOfCorrect()
{
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}

string DescriptionOfCorrect::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";
string DescriptionOfCorrect::isoAlphabet = "DECYHKR";   
float DescriptionOfCorrect::pKiso[7] = {-3.86,-4.25,-8.33,-10.0,6.0,10.5,12.4}; // Lehninger
float DescriptionOfCorrect::pKN = 9.69;
float DescriptionOfCorrect::pKC = 2.34;


void DescriptionOfCorrect::trainCorrect() {
  double piSum=0.0, dMSum=0.0;
  for(size_t ix=0; ix<psms.size(); ++ix) {
    piSum += psms[ix]->pI;
    dMSum += psms[ix]->massDiff;
  }
  if (psms.size()==0) {
    avgPI = 0.0; avgDM = 0.0;
  } else {
    avgPI = piSum/psms.size();
    avgDM = dMSum/psms.size();
  }
}
void DescriptionOfCorrect::setFeatures(PSMDescription* pPSM) {
  pPSM->features[DataSet::getDMFeatureNo()] = abs(pPSM->massDiff-avgDM);
  pPSM->features[DataSet::getPIFeatureNo()] = abs(pPSM->pI-avgPI);

}


double DescriptionOfCorrect::indexSum(const float* index, const string& peptide) {
  double sum = 0.0;
  string::const_iterator token = peptide.begin();
  for(;token != peptide.end();++token)
    sum += index[*token-'A'];
  return sum;
}

inline double DescriptionOfCorrect::indexN(const float *index, const string& peptide) {
  return index[peptide[0]-'A'];
}

inline double DescriptionOfCorrect::indexC(const float *index, const string& peptide) {
  return index[peptide[peptide.size()-1]-'A'];
}

inline double DescriptionOfCorrect::indexNC(const float *index, const string& peptide) {
  double n = max(0.0,indexN(index,peptide));
  double c = max(0.0,indexC(index,peptide));
  return n*c;
}

double* DescriptionOfCorrect::fillFeaturesIndex(const string& peptide, const float *index, double *features) {
  *(features++) = indexSum(index,peptide);
  *(features++) = indexN(index,peptide);
  *(features++) = indexC(index,peptide);
  *(features++) = indexNC(index,peptide);
  return features;
}

double* DescriptionOfCorrect::fillAAFeatures(const string& pep, double *feat) {
  // Overall amino acid composition features
  string::size_type pos = aaAlphabet.size();
  for (;pos--;) {feat[pos]=0.0;}
  for (string::const_iterator it=pep.begin();it!=pep.end();it++) {
    pos=aaAlphabet.find(*it);
    if (pos!=string::npos) feat[pos]++;
  }
  return feat+aaAlphabet.size();
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
    double NQ = 1/(1+pow(10,(pH-pKN))) - 1/(1+pow(10,(3.65-pKC)));
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

void DescriptionOfCorrect::fillFeaturesAllIndex(const string& peptide, double *features) {
  features = fillFeaturesIndex(peptide, krokhin_index, features);
  features = fillFeaturesIndex(peptide, hessa_index, features);
  features = fillFeaturesIndex(peptide, kytedoolittle_index, features);
  features = fillAAFeatures(peptide, features);
}


float DescriptionOfCorrect::krokhin_index['Z'-'A'+1] = 
         {0.8, 0.0, -0.8, -0.5, 0.0,  10.5, -0.9, -1.3, 8.4, 0.0, -1.9,9.6,5.8,
          -1.2,0.0,0.2,-0.9,-1.3,-0.8,0.4,0.0,5.0,11.0,0.0,4.0,0.0};
// negated hessa scale
float DescriptionOfCorrect::hessa_index['Z'-'A'+1] = 
         {-0.11,-0.0,0.13,-3.49,-2.68,0.32,-0.74,-2.06,0.60,0.0,-2.71,0.55,0.10,-2.05,
          0.0,-2.23,-2.36,-2.58,-0.84,-0.52,0.0,0.31,-0.30,0.0,-0.68,0.0};
float DescriptionOfCorrect::kytedoolittle_index['Z'-'A'+1] =
         {1.80,0.0,2.50,-3.50,-3.50,2.80,-0.40,-3.20,0.0,4.50,3.90,3.80,1.90,-3.50,
          0.0,-1.60,-3.50,-4.50,-0.80,-0.70,0.0,4.20,-0.90,0.0,-1.30,0.0};


