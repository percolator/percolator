/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SanityCheck.cpp,v 1.3 2008/05/23 22:12:20 lukall Exp $
 *******************************************************************************/
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"

// Class doing the sanity check for non SQT file condition
// In most sence a place holder as very little logic is build-in in this case

SanityCheck::SanityCheck() : initPositives(0),pTestset(NULL),pTrainset(NULL)
{
}

SanityCheck::~SanityCheck()
{
}

bool SanityCheck::overRule=false;
string SanityCheck::initWeightFN="";

int SanityCheck::getInitDirection(Scores * testset, Scores * trainset,Normalizer * pNorm,double *w,double test_fdr) {
  pTestset = testset; pTrainset = trainset; fdr = test_fdr;
  
  if (initWeightFN.size()>0) {
    C_DARRAY(ww,DataSet::getNumFeatures()+1)
    ifstream weightStream(initWeightFN.data(),ios::in);
    readWeights(weightStream,ww);
    weightStream.close();
    pNorm->normalizeweight(ww,w); 
    D_DARRAY(ww)
  } else {
    getDefaultDirection(w);
  }
  initPositives = pTestset->getInitDirection(fdr,w,false);
  return initPositives;
}

void SanityCheck::getDefaultDirection(double *w) {
  // Set init direction to be the most discriminative direction
  pTrainset->getInitDirection(fdr,w,true);
}


bool SanityCheck::validateDirection(double *w) {
  if (!pTestset) {
    cerr << "Wrongly set up of object SanityCheck" << endl;
    exit(-1);
  }
  int overFDR = pTestset -> calcScores(w,fdr);
  if (VERB>0) cerr << "Found " << overFDR << " peptides scoring over " << fdr*100 << "% FDR level on testset" << endl;
  if (overFDR<=0) {
     cerr << "No target score better than best decoy" << endl;
     resetDirection(w);
     return false;
  }
  if (initPositives>=overFDR) {
     cerr << "Less identifications after percolator processing than before processing" << endl;
     resetDirection(w);
     return false;
  }
  return true;
}

void SanityCheck::readWeights(istream & weightStream, double * w) {
  char buffer[1024],c;
  while (!(((c = weightStream.get())== '-') || (c >= '0' && c <= '9'))) {
    weightStream.getline(buffer,1024);
  }
  weightStream.getline(buffer,1024);
// Get second line containing raw features
  for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) {
    weightStream >> w[ix];
  }
}

void SanityCheck::resetDirection(double* w) {
  if (!overRule) {
    cerr << "Reseting score vector, using default vector" << endl;
    getDefaultDirection(w);  
  }
}


