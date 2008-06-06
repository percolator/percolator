/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.cpp,v 1.63 2008/06/06 17:13:32 lukall Exp $
 *******************************************************************************/
#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"
#include "PosteriorEstimator.h"
#include "ssl.h"

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score>other.score);}

inline bool operator<(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score<other.score);}

Scores::Scores()
{
    factor=1;
    neg=0;
    pos=0;
    posNow=0;
}

Scores::~Scores()
{
}

double Scores::pi0 = 0.9;

void Scores::merge(vector<Scores>& sv) {
  scores.clear();
  for(vector<Scores>::const_iterator a = sv.begin();a!=sv.end();a++)
    copy(a->begin(),a->end(),back_inserter(scores));
  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
}


void Scores::printRoc(string & fn){
 ofstream rocStream(fn.data(),ios::out);
 vector<ScoreHolder>::iterator it;
 for(it=scores.begin();it!=scores.end();it++) {
   rocStream << (it->label==-1?-1:1) << endl;
 }
 rocStream.close();
}	

double Scores::calcScore(const double * feat) const{
  register int ix=DataSet::getNumFeatures();
  register double score = w_vec[ix];
  for(;ix--;) {
  	score += feat[ix]*w_vec[ix];
  }
  return score;
}

ScoreHolder* Scores::getScoreHolder(const double *d){
  if (scoreMap.size()==0) {
    vector<ScoreHolder>::iterator it;
    for(it=scores.begin();it!=scores.end();it++) {
      scoreMap[it->featVec] = &(*it);
    }
  }
  return scoreMap[d];  
}


void Scores::fillFeatures(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio) {
  assert(ratio>0 && ratio < 1);
  int n = norm.getSize();
  int l = shuff.getSize();
  train.scores.clear();
  test.scores.clear();

  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  const double * featVec;
  while((featVec=shuffIter.getNext())!=NULL) {
    train.scores.push_back(ScoreHolder(.0,-1,featVec));
    test.scores.push_back(ScoreHolder(.0,-1,featVec));
  }
  while((featVec=normIter.getNext())!=NULL) {
    train.scores.push_back(ScoreHolder(.0,1,featVec));
    test.scores.push_back(ScoreHolder(.0,1,featVec));
  }
  train.pos=n;
  test.pos=n;
  train.neg=l;
  test.neg=l;
  train.factor = n/(double)l;
  test.factor = n/(double)l;
}

void Scores::fillFeatures(SetHandler& norm,SetHandler& shuff) {
  scores.clear();
  const double * featVec;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  while((featVec=normIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,1,featVec));
  }
  while((featVec=shuffIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,-1,featVec));
  }
  factor=norm.getSize()/(double)shuff.getSize();
}


void Scores::createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold) {
  train.resize(xval_fold);
  test.resize(xval_fold);
  for(unsigned int j=0;j<scores.size();j++) {
    for(unsigned int i=0;i<xval_fold;i++) {
      if(j%xval_fold==i) {
        test[i].scores.push_back(scores[j]);
      } else {
        train[i].scores.push_back(scores[j]);
      }
    }  
  }
  vector<ScoreHolder>::const_iterator it;
  for(unsigned int i=0;i<xval_fold;i++) {
  	train[i].pos=0;train[i].neg=0;
  	for(it=train[i].begin();it!=train[i].end();it++) {
      if (it->label==1) train[i].pos++;
      else train[i].neg++;
    }
    train[i].factor=train[i].pos/(double)train[i].neg;
    test[i].pos=0;test[i].neg=0;
  	for(it=test[i].begin();it!=test[i].end();it++) {
      if (it->label==1) test[i].pos++;
      else test[i].neg++;
  	}
    test[i].factor=test[i].pos/(double)test[i].neg;
  }
}

int Scores::calcScores(vector<double>& w,double fdr) {
  register unsigned int ix=DataSet::getNumFeatures()+1;
  if (w_vec.size()!=ix)
    w_vec.resize(ix);
  for(;ix--;) {
    w_vec[ix]=w[ix];
  }
  const double * features;
  vector<ScoreHolder>::iterator it = scores.begin();
  while(it!=scores.end()) {
    features = it->featVec;
  	it->score = calcScore(features);
    it++;
  }
  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  if (VERB>3) {
    cerr << "10 best scores and labels" << endl;
    unsigned int ix;
    for (ix=0;ix < 10;ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
    cerr << "10 worst scores and labels" << endl;
    for (ix=scores.size()-10;ix < scores.size();ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
  }
  int positives=0,nulls=0;
  double efp=0.0,q;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    if (q>pi0)
      q=pi0;
    it->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
  return posNow;
}

/*
double Scores::getQ(const double score) {
  unsigned int loIx = scores.size()-1, hiIx = 0;
  double lo = scores[loIx].score, hi = scores[hiIx].score;
  if (score<=lo) {return qVals[loIx];}
  if (score>=hi) {return qVals[hiIx];}
  int ix = 0, pos = -1;
  if (shortCut.empty()) {
    shortStep = (lo - hi)/(double)shortCutSize;
    while(++pos<shortCutSize) {
      double val = hi + shortStep*pos;
      while (scores[ix].score>=val)
        ++ix;
      shortCut.push_back(ix-1);
    }
    shortCut.push_back(loIx);   
  }
  pos = (int) ((score - hi)/shortStep);
  hiIx = shortCut[pos];
  if (shortCut[pos+1]==hiIx && hiIx<loIx) {
    loIx=hiIx+1;
  } else {
    loIx = shortCut[pos+1];
  }
  hi = scores[hiIx].score;
  lo = scores[loIx].score; 
  while((loIx-hiIx>1) && (qVals[loIx] > qVals[hiIx])) {
//    double d = (hiIx-loIx)/(hi-lo);
    ix = hiIx + (loIx-hiIx)/2;
    double sc = scores[ix].score;
    if (sc<score) {
      lo = sc;
      loIx = ix;
    } else {
      hi = sc;
      hiIx = ix;    
    }
  }
  return qVals[loIx];
}
*/


void Scores::generateNegativeTrainingSet(AlgIn& data,const double cneg) {
  unsigned int ix1=0,ix2=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==-1) {
      data.vals[ix2]=scores[ix1].featVec;
      data.Y[ix2]=-1;
      data.C[ix2++]=cneg;
    }
  }
  data.negatives=ix2;
}


void Scores::generatePositiveTrainingSet(AlgIn& data,const double fdr,const double cpos) {
  unsigned int ix1=0,ix2=data.negatives,p=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==1) {
      if (fdr<scores[ix1].q) {
        posNow=p;
        break;
      }
      data.vals[ix2]=scores[ix1].featVec;
      data.Y[ix2]=1;
      data.C[ix2++]=cpos;
      ++p;
    }
  }
  data.positives=p;
  data.m=ix2;
}

int Scores::getInitDirection(const double fdr, vector<double>& direction, bool findDirection) {
  int bestPositives = -1;
  int bestFeature =-1;
  bool lowBest = false;
  
  if (findDirection) { 
    for (int featNo=0;featNo<DataSet::getNumFeatures();featNo++) {
      vector<ScoreHolder>::iterator it = scores.begin();
      while(it!=scores.end()) {
        it->score = it->featVec[featNo];
        it++;
      }
      sort(scores.begin(),scores.end());
      for (int i=0;i<2;i++) {
        int positives=0,nulls=0;
        double efp=0.0,q;
        for(it=scores.begin();it!=scores.end();it++) {
          if (it->label!=-1)
            positives++;
          if (it->label==-1) {
            nulls++;
            efp=pi0*nulls*factor;
          }
          if (positives)
            q=efp/(double)positives;
          else
            q=pi0;
          if (fdr<=q) {
            if (positives>bestPositives && scores.begin()->score!=it->score) {
              bestPositives=positives;
              bestFeature = featNo;
              lowBest = (i==0);
            }
            if (i==0) {
              reverse(scores.begin(),scores.end());
            }
            break;
          }
        }
      }
    }
    for (int ix=DataSet::getNumFeatures();ix--;) {
      direction[ix]=0;
    }
    direction[bestFeature]=(lowBest?-1:1);
    if (VERB>1) {
      cerr << "Selected feature number " << bestFeature +1 << " as initial search direction, could separate " << 
              bestPositives << " positives in that direction" << endl;
    }
  } else {
    bestPositives = calcScores(direction,fdr);
    if (VERB>1) {
      cerr << "Found " << 
              bestPositives << " positives in the initial search direction" << endl;
    }    
  }
  return bestPositives;
}

double Scores::estimatePi0() {
  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));

  // Estimate pi0
  pi0 = PosteriorEstimator::estimatePi0(combined);
  return pi0;

}

void Scores::calcPep() {

  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));
  vector<double> peps;                                                                                                                  
  
  // Logistic regression on the data
  PosteriorEstimator::estimatePEP(combined,pi0,peps);

  size_t pix=0;
  for(size_t ix=0; ix<scores.size(); ix++) {
    scores[ix].pep = peps[pix];
    if(scores[ix].label==1)
      ++pix;
  }
}
