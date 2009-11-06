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
#include "svm.h"
#include "Globals.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "DescriptionOfCorrect.h"
#include "MassHandler.h"
#include "PSMDescription.h"

static double REGRESSION_C = 5.0;

DescriptionOfCorrect::DescriptionOfCorrect()
{
  numRTFeat = 0;
  model = NULL;
  c=REGRESSION_C;
  gamma = 50;
  epsilon = 1e-3;
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}


string DescriptionOfCorrect::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";
string DescriptionOfCorrect::isoAlphabet = "DECYHKR";
string DescriptionOfCorrect::ptmAlphabet = "#*@";
float DescriptionOfCorrect::pKiso[7] = {-3.86,-4.25,-8.33,-10.0,6.0,10.5,12.4}; // Lehninger
float DescriptionOfCorrect::pKN = 9.69;
float DescriptionOfCorrect::pKC = 2.34;
bool DescriptionOfCorrect::doKlammer=false;
unsigned int DescriptionOfCorrect::docFeatures = 15;

void DescriptionOfCorrect::calcRegressionFeature(PSMDescription &psm) {
  string peptide = psm.getPeptide();
  string::size_type pos1 = peptide.find('.');
  string::size_type pos2 = peptide.find('.',++pos1);
  string pep = peptide.substr(pos1,pos2-pos1);
  psm.pI = isoElectricPoint(pep);
  if (psm.getRetentionFeatures()) {
    fillFeaturesAllIndex(pep, psm.getRetentionFeatures());
  }
//  cout <<  peptide << " " << pep << " " << retentionFeatures[0] << endl;
}

double DescriptionOfCorrect::testRetention(vector<PSMDescription>& testset) {
  double rms=0.0;
  for(size_t ix1=0;ix1<testset.size();ix1++) {
    double diff = estimateRT(testset[ix1].retentionFeatures)-testset[ix1].retentionTime;
	rms += diff*diff;
  }
  return rms/testset.size();
}

void DescriptionOfCorrect::trainRetention(vector<PSMDescription>& trainset, const double C, const double gamma, const double epsilon) {
  if (psms.size()>500) {
    numRTFeat = totalNumRTFeatures();
  } else {
    numRTFeat = minimumNumRTFeatures();
  }
  svm_parameter param;
  param.svm_type = EPSILON_SVR;
  param.kernel_type = RBF;
  param.degree = 3;
  param.gamma = 1/(double)psms.size()*gamma;
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = C;
  param.eps = epsilon; //1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;

  svm_problem data;
  data.l=trainset.size();
  data.x=new svm_node[data.l];
  data.y=new double[data.l];

  for(size_t ix1=0;ix1<trainset.size();ix1++) {
    data.x[ix1].values=trainset[ix1].retentionFeatures;
    data.x[ix1].dim=numRTFeat;
    data.y[ix1]=trainset[ix1].retentionTime;
  }
  svm_model* m = svm_train(&data, &param);
  copyModel(m);
  delete[] data.x;
  delete[] data.y;
  svm_destroy_model(m);
}

void DescriptionOfCorrect::trainCorrect() {
  // Get rid of redundant peptides
  sort(psms.begin(),psms.end(),less<PSMDescription>());
  psms.resize(distance(psms.begin(),unique(psms.begin(),psms.end())));

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

  // Train teyrntion time regressor
  size_t test_frac = 4u;
  if (psms.size()>test_frac*10u) {
	// If we got enough data, calibrate gamma and C by leaving out a testset
    vector<PSMDescription> train,test;
    for(size_t ix=0; ix<psms.size(); ++ix) {
      if (ix%test_frac==0) {
    	  test.push_back(psms[ix]);
      } else {
    	  train.push_back(psms[ix]);
      }
    }
    double bestRms = 1e100;
    double gammaV[3] = {gamma/2,gamma,gamma*2};
    double cV[3] = {c/2,c,c*2};
    double epsilonV[3] = {epsilon/2,epsilon,epsilon*2};
    for (double* gammaNow=&gammaV[0];gammaNow!=&gammaV[3];gammaNow++){
        for (double* cNow=&cV[0];cNow!=&cV[3];cNow++){
            for (double* epsilonNow=&epsilonV[0];epsilonNow!=&epsilonV[3];epsilonNow++){
        	  trainRetention(train,*cNow,(*gammaNow),*epsilonNow);
        	  double rms=testRetention(test);
        	  if (rms<bestRms) {
        		  c=*cNow;gamma=*gammaNow;epsilon=*epsilonNow;
        		  bestRms=rms;
        	  }
            }
        }
    }
    // cerr << "CV selected gamma=" << gamma << " and C=" << c << endl;
  }
  trainRetention(psms,c,gamma,epsilon);
  if (VERB>2) cerr << "Description of correct recalibrated, avg pI=" << avgPI << " avg dM=" << avgDM << endl;

}
void DescriptionOfCorrect::setFeatures(PSMDescription& psm) {
  assert(DataSet::getFeatureNames().getDocFeatNum()>0);
  psm.predictedTime=estimateRT(psm.retentionFeatures);
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

void DescriptionOfCorrect::copyModel(svm_model* from) {
  if (model!=NULL)
    svm_destroy_model(model);
  model = (svm_model *) malloc(sizeof(svm_model));
  model->param=from->param;
  model->nr_class = from->nr_class;
  model->l = from->l;			// total #SV
// _DENSE_REP
  model->SV = (svm_node *)malloc(sizeof(svm_node)*model->l);
  for(int i=0;i<model->l;++i) {
    model->SV[i].dim = from->SV[i].dim;
    model->SV[i].values = (double *)malloc(sizeof(double)*model->SV[i].dim);
    memcpy(model->SV[i].values,from->SV[i].values,sizeof(double)*model->SV[i].dim);
  }
  model->sv_coef = (double **) malloc(sizeof(double *)*(model->nr_class-1));
  for(int i=0;i<model->nr_class-1;++i) {
    model->sv_coef[i] = (double *) malloc(sizeof(double)*(model->l));
    memcpy(model->sv_coef[i],from->sv_coef[i],sizeof(double)*(model->l));
  }
  int n = model->nr_class * (model->nr_class-1)/2;
  if (from->rho) {
    model->rho = (double *) malloc(sizeof(double)*n);
    memcpy(model->rho,from->rho,sizeof(double)*n);
  } else {
    model->rho=NULL;
  }
  if (from->probA) {
    model->probA = (double *) malloc(sizeof(double)*n);
    memcpy(model->probA,from->probA,sizeof(double)*n);
  } else { model->probA=NULL; }
  if (from->probB) {
    model->probB = (double *) malloc(sizeof(double)*n);
    memcpy(model->probB,from->probB,sizeof(double)*n);
  } else { model->probB=NULL; }

  assert(from->label==NULL);
  assert(from->nSV==NULL);

  // for classification only
  model->label=NULL;		// label of each class (label[k])
  model->nSV=NULL;		// number of SVs for each class (nSV[k])
  // nSV[0] + nSV[1] + ... + nSV[k-1] = l
  // XXX
  model->free_sv=1;         // 1 if svm_model is created by svm_load_mode
	                        // 0 if svm_model is created by svm_train
}

double DescriptionOfCorrect::estimateRT(double * features) {
  svm_node node;
  node.values = features;
  node.dim = numRTFeat;
  return svm_predict(model,&node);

/*
  register int ix = rtW.size()-1;
  double sum = rtW[ix];
  for(;ix--;)
    sum += rtW[ix]*features[ix];
  return sum;
  */
}


double DescriptionOfCorrect::indexSum(const float* index, const string& peptide) {
  double sum = 0.0;
  string::const_iterator token = peptide.begin();
  for(;token != peptide.end();++token)
    sum += index[*token-'A'];
  return sum;
}

double DescriptionOfCorrect::indexAvg(const float* index, const string& peptide) {
  double sum = 0.0;
  string::const_iterator token = peptide.begin();
  for(;token != peptide.end();++token)
    sum += index[*token-'A'];
  return sum/(double)peptide.size();
}

double DescriptionOfCorrect::indexNearestNeigbourPos(const float* index, const string& peptide) {
  double sum = 0.0;
  for(unsigned int ix = 0; ix < peptide.size();++ix) {
    if (peptide[ix]=='R' || peptide[ix]=='K') {
      if (ix>0)
        sum += max(0.0f,index[peptide[ix-1]-'A']);
      if (ix<peptide.size()-1)
        sum += max(0.0f,index[peptide[ix+1]-'A']);
    }
  }
  return sum;
}

double DescriptionOfCorrect::indexNearestNeigbourNeg(const float* index, const string& peptide) {
  double sum = 0.0;
  for(unsigned int ix = 0; ix < peptide.size();++ix) {
    if (peptide[ix]=='D' || peptide[ix]=='E') {
      if (ix>0)
        sum += max(0.0f,index[peptide[ix-1]-'A']);
      if (ix<peptide.size()-1)
        sum += max(0.0f,index[peptide[ix+1]-'A']);
    }
  }
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

double* DescriptionOfCorrect::indexPartialSum(const float* index, const string& peptide, const size_t win, double *features) {
  double sum = 0.0;
  size_t window = min(win,peptide.size()-1);
  string::const_iterator lead = peptide.begin(), lag = peptide.begin();
  for(;lead != (peptide.begin() + window);++lead)
    sum += index[*lead-'A'];
  double minS = sum, maxS = sum;
  for(;lead != peptide.end();++lead,++lag) {
    sum += index[*lead-'A'];
    sum -= index[*lag-'A'];
    minS = min(sum,minS);
    maxS = max(sum,maxS);
  }
  *(features++) = maxS;
  *(features++) = minS;
  return features;
}


double* DescriptionOfCorrect::fillFeaturesIndex(const string& peptide, const float *index, double *features) {
  *(features++) = indexSum(index,peptide);
  *(features++) = indexAvg(index,peptide);
  *(features++) = indexN(index,peptide);
  *(features++) = indexC(index,peptide);
  *(features++) = indexNearestNeigbourPos(index,peptide);
  *(features++) = indexNearestNeigbourNeg(index,peptide);
  features = indexPartialSum(index,peptide,3,features);
  features = indexPartialSum(index,peptide,5,features);
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
          cerr << psms[i].retentionFeatures[j] << "\t";
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
    *(features++) = DataSet::isEnz(Ct,'A');
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
