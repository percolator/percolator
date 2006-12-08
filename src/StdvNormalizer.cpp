#include <vector>
#include <iostream>
#include <math.h>
#include <set>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "SetHandler.h"
#include "Globals.h"

StdvNormalizer::StdvNormalizer()
{
  avg.resize(DataSet::getNumFeatures(),0.0);
  stdv.resize(DataSet::getNumFeatures(),0.0);
}

StdvNormalizer::~StdvNormalizer()
{
}

void StdvNormalizer::normalize(const double *in,double* out){
  for (int ix=0;ix<DataSet::getNumFeatures();ix++) {
  	out[ix]=(in[ix]-avg[ix])/stdv[ix];
  }
}

void StdvNormalizer::unnormalizeweight(const double *in,double* out){
  double sum = 0;
  int i=0;
  for (;i<DataSet::getNumFeatures();i++) {
  	out[i]=in[i]/stdv[i];
  	sum+=avg[i]*in[i]/stdv[i];
  }
  out[i]=in[i]-sum;
}

void StdvNormalizer::normalizeweight(const double *in,double* out){
  double sum = 0;
  int i=0;
  for (;i<DataSet::getNumFeatures();i++) {
  	out[i]=in[i]*stdv[i];
  	sum+=avg[i]*in[i];
  }
  out[i]=in[i]+sum;
}

void StdvNormalizer::setSet(set<DataSet *> & setVec){
  double n=0.0;
  const double * features;
  int ix;
  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
    avg[ix]=0.0;
    stdv[ix]=0.0;
  }
  set<DataSet *>::iterator it;
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((features=(*it)->getNext(ixPos))!=NULL) {
	  n++;
	  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
	    avg[ix]+=features[ix];
	  }
    }
  }
  if (VERB>2) { 
    cerr.precision(2);
    cerr << "Normalization factors" << endl
    << "Type\t" << DataSet::getFeatureNames() << endl << "Avg ";
  }
  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
  	if (n>0.0)
     avg[ix]/=n;
     if (VERB>2) cerr << "\t" << avg[ix]; 
  }
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((features=(*it)->getNext(ixPos))!=NULL) {
      for (ix=0;ix<DataSet::getNumFeatures();ix++) {
        if (!isfinite(features[ix]))
          cerr << "Reached strange feature with val=" << features[ix] << " at row=" << ix << ", col=" << ixPos << endl;
        double d = features[ix]-avg[ix];
        stdv[ix]+=d*d;
      }
    }
  }
  if (VERB>2) cerr << endl << "Stdv"; 
  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
    if (stdv[ix]<=0 || n==0) {
      stdv[ix]=1.0;
    } else {
  	  stdv[ix]=sqrt(stdv[ix]/n);
  	  if (VERB>2) cerr << "\t" << stdv[ix]; 
    }
  }
  if (VERB>2) cerr << endl; 
}
