/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: UniNormalizer.cpp,v 1.17 2009/01/04 22:49:31 lukall Exp $
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
}

UniNormalizer::~UniNormalizer()
{
}

void UniNormalizer::unnormalizeweight(const vector<double>& in,vector<double>& out){
  double sum = 0;
  unsigned int i=0;
  for (;i<numFeatures;i++) {
  	out[i]=in[i]/div[i];
  	sum += sub[i]*in[i]/div[i];
  }
  out[i]=in[i]-sum;
}

void UniNormalizer::normalizeweight(const vector<double>& in, vector<double>& out){
  double sum = 0;
  size_t i=0;
  for (;i<numFeatures;i++) {
  	out[i]=in[i]*div[i];
  	sum+=sub[i]*in[i];
  }
  out[i]=in[i]+sum;
}

void UniNormalizer::setSet(set<DataSet *> &setVec, size_t nf, size_t nrf){
  numFeatures = nf; numRetentionFeatures=nrf;
  sub.resize(nf+nrf,0.0); div.resize(nf+nrf,0.0);
  vector<double> mins(nf+nrf,1e+100), maxs(nf+nrf,-1e+100);

  double * features;
  PSMDescription* pPSM;
  size_t ix;
  set<DataSet *>::iterator it;
  for (it=setVec.begin();it!=setVec.end();++it) {
    int ixPos=-1;
    while((pPSM = (*it)->getNext(ixPos))!=NULL) {
      features = pPSM->features;
	  for (ix=0;ix<numFeatures;ix++) {
	    mins[ix]=min(features[ix],mins[ix]);
	    maxs[ix]=max(features[ix],maxs[ix]);
      }
      features = pPSM->retentionFeatures;
      for (;ix<numFeatures+numRetentionFeatures;++ix) {
        mins[ix]=min(features[ix-numFeatures],mins[ix]);
        maxs[ix]=max(features[ix-numFeatures],maxs[ix]);
      }
	}
  }
  for (ix=0;ix<numFeatures+numRetentionFeatures;++ix) {
  	sub[ix]=mins[ix];
  	div[ix]=maxs[ix]-mins[ix];
  	if (div[ix]<=0)
  	  div[ix]=1.0;
  }
}
