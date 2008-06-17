/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: StdvNormalizer.cpp,v 1.20 2008/06/17 23:21:44 lukall Exp $
 *******************************************************************************/
#include <vector>
#include <iostream>
#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif
#include <math.h>
#include <set>
#include <vector>
#include <string>
using namespace std;
#include "PSMDescription.h"
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

void StdvNormalizer::unnormalizeweight(const vector<double>& in, vector<double>& out){
  double sum = 0;
  int i=0;
  for (;i<DataSet::getNumFeatures();i++) {
  	out[i]=in[i]/stdv[i];
  	sum+=avg[i]*in[i]/stdv[i];
  }
  out[i]=in[i]-sum;
}

void StdvNormalizer::normalizeweight(const vector<double>& in, vector<double>& out){
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
  double * features;
  PSMDescription* pPSM;
  int ix;
  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
    avg[ix]=0.0;
    stdv[ix]=0.0;
  }
  set<DataSet *>::iterator it;
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((pPSM=(*it)->getNext(ixPos))!=NULL) {
      features = pPSM->features;
	  n++;
	  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
	    avg[ix]+=features[ix];
	  }
    }
  }
  if (VERB>2) { 
    cerr.precision(2);
    cerr << "Normalization factors" << endl
    << "Type\t" << DataSet::getFeatureNames().getFeatureNames() << endl << "Avg ";
  }
  for (ix=0;ix<DataSet::getNumFeatures();ix++) {
  	if (n>0.0)
     avg[ix]/=n;
     if (VERB>2) cerr << "\t" << avg[ix]; 
  }
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((pPSM=(*it)->getNext(ixPos))!=NULL) {
      features = pPSM->features;
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
    }
    if (VERB>2) cerr << "\t" << stdv[ix]; 
  }
  if (VERB>2) cerr << endl; 
}
