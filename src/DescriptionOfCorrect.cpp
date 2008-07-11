#include <vector>
#include <cmath>
#include <assert.h>
#include "ssl.h"
#include "Globals.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "DescriptionOfCorrect.h"

DescriptionOfCorrect::DescriptionOfCorrect()
{
  numRTFeat = 0;
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}

static double REGRESSION_C = 0.1; 

string DescriptionOfCorrect::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";
string DescriptionOfCorrect::isoAlphabet = "DECYHKR";   
float DescriptionOfCorrect::pKiso[7] = {-3.86,-4.25,-8.33,-10.0,6.0,10.5,12.4}; // Lehninger
float DescriptionOfCorrect::pKN = 9.69;
float DescriptionOfCorrect::pKC = 2.34;
bool DescriptionOfCorrect::doIsotopeMass=false;


void DescriptionOfCorrect::trainRetention() {
  if (psms.size()>500) {
    numRTFeat = totalNumRTFeatures();
  } else {
    numRTFeat = 12;  
  }
  AlgIn data(psms.size(),numRTFeat);
  data.m = psms.size();
  unsigned int ix1=0;
  for(ix1=0;ix1<psms.size();ix1++) {
    data.vals[ix1]=psms[ix1]->retentionFeatures;
    data.Y[ix1]=psms[ix1]->retentionTime;
    data.C[ix1]=REGRESSION_C;
  }
  struct options *Options = new options;
  Options->lambda=1.0;
  Options->lambda_u=1.0;
  Options->epsilon=EPSILON;
  Options->cgitermax=CGITERMAX;
  Options->mfnitermax=MFNITERMAX;
    
  struct vector_double *Weights = new vector_double;
  Weights->d = numRTFeat+1;
  Weights->vec = new double[Weights->d];
//    for(int ix=0;ix<Weights->d;ix++) Weights->vec[ix]=w[ix];
  for(int ix=0;ix<Weights->d;ix++) Weights->vec[ix]=0;
  rtW.resize(numRTFeat+1);
    
  struct vector_double *Outputs = new vector_double;
  Outputs->vec = new double[psms.size()];
  Outputs->d = psms.size();
  for(int ix=0;ix<Outputs->d;ix++) Outputs->vec[ix]=0;

  L2_SVM_MFN(data,Options,Weights,Outputs);

  for(int i= rtW.size();i--;)
    rtW[i]=Weights->vec[i];

  delete [] Weights->vec;
  delete Weights;
  delete [] Outputs->vec;
  delete Outputs;
  delete Options;
}

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
  // print_10features();
  trainRetention();
  if (VERB>2) cerr << "Description of correct recalibrated, avg pI=" << avgPI << " amd avg dM=" << avgDM << endl;  
  
}
void DescriptionOfCorrect::setFeatures(PSMDescription* pPSM) {
  assert(DataSet::getFeatureNames().getDocFeatNum()>0);
  size_t docFeatNum = DataSet::getFeatureNames().getDocFeatNum();
  pPSM->features[docFeatNum] = Normalizer::getNormalizer()->normalize(abs(pPSM->pI-avgPI),docFeatNum);
  pPSM->features[docFeatNum+1] = Normalizer::getNormalizer()->normalize(deltadeltaMass(pPSM->massDiff),docFeatNum+1);
  pPSM->features[docFeatNum+2] = Normalizer::getNormalizer()->normalize(abs(pPSM->retentionTime-estimateRT(pPSM->retentionFeatures)),docFeatNum+2);
//  pPSM->features[DataSet::getFeatureNames().getDocFeatNum()+2] = 0.0;

}

inline double DescriptionOfCorrect::estimateRT(double * features) {
  register int ix = numRTFeat;
  double sum = rtW[ix];
  for(;ix--;) 
    sum += rtW[ix]*features[ix];
  return sum;
}

inline double DescriptionOfCorrect::deltadeltaMass(double dm) {
  double ddm = dm - avgDM;
  if (!doIsotopeMass) 
    return abs(ddm);
  double isoddm = abs(ddm-1);
  for(int isotope=0;isotope<5;++isotope) {
    isoddm = min(isoddm,abs(ddm+isotope));
  }  
  return isoddm;
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

void DescriptionOfCorrect::print_10features() {
   for(int i=0;i<10;i++) {
       for(size_t j=0;j<totalNumRTFeatures();j++) {
          cerr << psms[i]->retentionFeatures[j] << "\t";
       }
       cerr << endl;
   }
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


