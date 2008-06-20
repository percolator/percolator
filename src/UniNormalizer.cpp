/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: UniNormalizer.cpp,v 1.15 2008/06/20 23:55:35 lukall Exp $
 *******************************************************************************/
#include <set>
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "UniNormalizer.h"
#include "SetHandler.h"

UniNormalizer::UniNormalizer()
{
  sub.resize(FeatureNames::getNumFeatures(),0.0);
  div.resize(FeatureNames::getNumFeatures(),0.0);
}

UniNormalizer::~UniNormalizer()
{
}

void UniNormalizer::normalize(const double *in,double* out){
  for (unsigned int ix=0;ix<FeatureNames::getNumFeatures();ix++) {
  	out[ix]=(in[ix]-sub[ix])/div[ix];
  }
}

void UniNormalizer::unnormalizeweight(const vector<double>& in,vector<double>& out){
  double sum = 0;
  unsigned int i=0;
  for (;i<FeatureNames::getNumFeatures();i++) {
  	out[i]=in[i]/div[i];
  	sum += sub[i]*in[i]/div[i];
  }
  out[i]=in[i]-sum;
}

void UniNormalizer::normalizeweight(const vector<double>& in, vector<double>& out){
  double sum = 0;
  unsigned int i=0;
  for (;i<FeatureNames::getNumFeatures();i++) {
  	out[i]=in[i]*div[i];
  	sum+=sub[i]*in[i];
  }
  out[i]=in[i]+sum;
}

void UniNormalizer::setSet(set<DataSet *> &setVec){
  double * features;
  PSMDescription* pPSM;
  unsigned int ix;
  vector<double> mins(FeatureNames::getNumFeatures(),1e+100);
  vector<double> maxs(FeatureNames::getNumFeatures(),-1e+100);
  set<DataSet *>::iterator it;
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((pPSM = (*it)->getNext(ixPos))!=NULL) {
      features = pPSM->features;
	  for (ix=0;ix<FeatureNames::getNumFeatures();ix++) {
	    mins[ix]=min(features[ix],mins[ix]);
	    maxs[ix]=max(features[ix],maxs[ix]);
      }
	}
  }
  for (ix=0;ix<FeatureNames::getNumFeatures();ix++) {
  	sub[ix]=mins[ix];
  	div[ix]=maxs[ix]-mins[ix];
  	if (div[ix]<=0)
  	  div[ix]=1.0;
  }
}
