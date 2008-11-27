#include <vector>
#include <cstring>
#include <cmath>
#include <assert.h>
#include "svm.h"
#include "Globals.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "DescriptionOfCorrect.h"

DescriptionOfCorrect::DescriptionOfCorrect()
{
  numRTFeat = 0;
  model = NULL;
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}

static double REGRESSION_C = 1.0; 

string DescriptionOfCorrect::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";
string DescriptionOfCorrect::isoAlphabet = "DECYHKR";   
string DescriptionOfCorrect::ptmAlphabet = "#*@";   
float DescriptionOfCorrect::pKiso[7] = {-3.86,-4.25,-8.33,-10.0,6.0,10.5,12.4}; // Lehninger
float DescriptionOfCorrect::pKN = 9.69;
float DescriptionOfCorrect::pKC = 2.34;
bool DescriptionOfCorrect::doIsotopeMass=false;
bool DescriptionOfCorrect::doKlammer=false;

void DescriptionOfCorrect::trainRetention() {
  if (psms.size()>500) {
    numRTFeat = totalNumRTFeatures();
  } else {
    numRTFeat = minimumNumRTFeatures();  
  }
  svm_parameter param;  
  param.svm_type = EPSILON_SVR;
  param.kernel_type = RBF;
  param.degree = 3;
  param.gamma = 1/(double)psms.size();
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = REGRESSION_C;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  
  svm_problem data; 
  data.l=psms.size();
  data.x=new svm_node[data.l];
  data.y=new double[data.l];

  for(size_t ix1=0;ix1<psms.size();ix1++) {
    data.x[ix1].values=psms[ix1]->retentionFeatures;
    data.x[ix1].dim=numRTFeat;
    data.y[ix1]=psms[ix1]->retentionTime;
  }
  svm_model* m = svm_train(&data, &param);
  copyModel(m);
  delete[] data.x;
  delete[] data.y;
  svm_destroy_model(m);
}
/*
void DescriptionOfCorrect::trainRetention() {
  if (psms.size()>500) {
    numRTFeat = totalNumRTFeatures();
  } else {
    numRTFeat = minimumNumRTFeatures();  
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
*/

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

double DescriptionOfCorrect::indexAvg(const float* index, const string& peptide) {
  double sum = 0.0;
  string::const_iterator token = peptide.begin();
  for(;token != peptide.end();++token)
    sum += index[*token-'A'];
  return sum/(double)peptide.size();
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

double* DescriptionOfCorrect::indexPartialSum(const float* index, const string& peptide, const size_t window, double *features) {
  double sum = 0.0;
  string::const_iterator lead = peptide.begin(), lag = peptide.begin();
  for(;lead != (peptide.begin() + window);++lead)
    sum += index[*lead-'A'];
  double minS = sum, maxS = sum;
  for(;lead != peptide.end();++lead,lag++) {
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
//  *(features++) = indexNC(index,peptide);
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

void DescriptionOfCorrect::fillFeaturesAllIndex(const string& pep, double *features) {
  unsigned int ptms=0;
  string peptide = pep;
  string::size_type posP = 0, pos = 0;
  while (posP < ptmAlphabet.length()) {
    while ( (pos = peptide.find(ptmAlphabet[posP], pos)) != string::npos ) {
      peptide.replace( pos, 1, "");
      ++pos;
      ++ptms;
    }
    ++posP;
  }  

  if(!doKlammer) {
    *(features++) = (double) ptms;    
    features = fillFeaturesIndex(peptide, krokhin_index, features);
    features = fillFeaturesIndex(peptide, hessa_index, features);
    features = fillFeaturesIndex(peptide, kytedoolittle_index, features);
    *(features++) = peptide.size();    
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
float DescriptionOfCorrect::aa_weights['Z'-'A'+1] = {
  71.0788,0, 103.1448,115.0886,129.1155,147.1766,57.052,137.1412,113.1595,0,128.1742,113.1595,131.1986,114.1039,0,97.1167,128.1308,156.1876,87.0782,101.1051,0,99.1326,186.2133,0,163.1760,0};
