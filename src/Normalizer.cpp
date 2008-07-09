/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Normalizer.cpp,v 1.18 2008/07/09 00:54:19 lukall Exp $
 *******************************************************************************/
#include <assert.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "UniNormalizer.h"
#include "Globals.h"

Normalizer::Normalizer()
{
}

Normalizer::~Normalizer()
{
}

void Normalizer::normalizeSet(set<DataSet *> & setVec) {
  double * features;
  PSMDescription* pPSM;
  set<DataSet *>::iterator it;
  if (VERB>4) {
    cerr << "First 10 feature vectors before normalization" << endl;
    cerr.precision(3);
    it=setVec.begin();
    int ixPos=-1;
    cerr << "Label of this set is " << (*it)->getLabel() << endl;
    while((pPSM=(*it)->getNext(ixPos))!=NULL && ixPos < 10) {
      features = pPSM->features;
      for (unsigned int a=0; a<numFeatures; ++a) {
        cerr << features[a] << " ";
      }
      cerr << endl;
    }
    
  }
  for (it=setVec.begin();it!=setVec.end();++it) {
      int ixPos=-1;
      while((pPSM=(*it)->getNext(ixPos))!=NULL) {
        features = pPSM->features;
        normalize(features,features,0,numFeatures);
        features = pPSM->retentionFeatures;
        normalize(features,features,numFeatures,numRetentionFeatures);
      }
  }
  if (VERB>4) {
    cerr << "First 10 feature vectors after normalization" << endl;
    cerr.precision(3);
    it=setVec.begin();
    int ixPos=-1;
    while((pPSM=(*it)->getNext(ixPos))!=NULL && ixPos < 10) {
      features = pPSM->features;     
      for (unsigned int a=0;a<numFeatures;++a) {
        cerr << features[a] << " ";
      }
      cerr << endl;
    }
    
  }
}

void Normalizer::normalize(const double *in,double* out, size_t offset, size_t numFeatures){
  for (unsigned int ix=0;ix<numFeatures;++ix) {
    out[ix]=(in[ix]-sub[offset+ix])/div[offset+ix];
  }
}

int Normalizer::subclass_type = STDV;

Normalizer * Normalizer::getNew(){
	if (subclass_type == UNI)
	  return new UniNormalizer();
	return new StdvNormalizer();
}

void Normalizer::setType(int type){
	assert(type==UNI || type==STDV);
	subclass_type = type;
}
