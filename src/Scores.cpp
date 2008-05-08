/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.cpp,v 1.57 2008/05/08 00:22:35 lukall Exp $
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
    shortStep=0;
    w_vec = NULL;
}

Scores::~Scores()
{
  if (w_vec)
    delete w_vec;
  w_vec = NULL;
}

double Scores::pi0 = 0.9;

void Scores::printRoc(string & fn){
 ofstream rocStream(fn.data(),ios::out);
 vector<ScoreHolder>::iterator it;
 for(it=scores.begin();it!=scores.end();it++) {
   rocStream << (it->label==-1?-1:1) << endl;
 }
 rocStream.close();
}	

double Scores::calcScore(const double *feat) const{
  register int ix=DataSet::getNumFeatures();
  register double score = w_vec[ix];
  for(;ix--;) {
  	score += feat[ix]*w_vec[ix];
  }
  return score;
}

void Scores::fillFeatures(Scores& train,Scores& thresh,Scores& test,SetHandler& norm,SetHandler& shuff,
                             const double trainRatio,const double testRatio) {
  assert(trainRatio >0 && testRatio >0 && trainRatio+testRatio < 1);
  int n = norm.getSize();
  int k = (int)(shuff.getSize()*trainRatio);
  int m = (int)(shuff.getSize()*testRatio);
  int l = shuff.getSize() - k - m;
  ScoreHolder s;
  train.scores.resize(n+k,s);
  train.qVals.resize(n+k,-1e200); 
  test.scores.resize(n+m,s);
  test.qVals.resize(n+m,-1e200); 
  thresh.scores.resize(n+l,s);
  thresh.qVals.resize(n+l,-1e200); 
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  int ix1=0,ix2=0,ix3=0;
  const double * featVec;
  while((featVec=shuffIter.getNext())!=NULL) {
    int remain = k+m+l - ix1 - ix2 -ix3;
    int sel = (rand() % remain);
    if (sel<(k-ix1)) {
//    if (((int)(ix1+ix2+ix3+1)*trainRatio)>=ix1+1) {
      train.scores[ix1].label=-1;
      train.scores[ix1].featVec=featVec;
      ++ix1;
    } else if (sel<(k-ix1)+(m-ix3)) {
//    } else if (((int)(ix1+ix2+ix3+1)*testRatio)>=ix3+1) {
      test.scores[ix3].label=-1;
      test.scores[ix3].featVec=featVec;
      ++ix3;    
    } else {
      thresh.scores[ix2].label=-1;
      thresh.scores[ix2].featVec=featVec;
      ++ix2;    
    }
  }
  assert(ix1==k);
  assert(ix2==l);
  assert(ix3==m);
  while((featVec=normIter.getNext())!=NULL) {
    train.scores[ix1].label=1;
    train.scores[ix1].featVec=featVec;
    ++ix1;
    thresh.scores[ix2].label=1;
    thresh.scores[ix2].featVec=featVec;
    ++ix2;
    test.scores[ix3].label=1;
    test.scores[ix3].featVec=featVec;
    ++ix3;
  }
  assert(ix1==n+k);
  assert(ix2==n+l);
  assert(ix3==n+m);
  train.pos=n;
  thresh.pos=n;
  test.pos=n;
  train.neg=k;
  thresh.neg=l;
  test.neg=m;
  train.factor = n/(double)k;
  thresh.factor = n/(double)l;
  test.factor = n/(double)m;
}

void Scores::fillFeatures(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio) {
  assert(ratio>0 && ratio < 1);
  int n = norm.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  ScoreHolder s;
  train.scores.resize(n+k,s);
  train.qVals.resize(n+k,-1e200); 
  test.scores.resize(n+l,s);
  test.qVals.resize(n+l,-1e200); 

  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  int ix1=0,ix2=0;
  const double * featVec;
  while((featVec=shuffIter.getNext())!=NULL) {
    int remain = k + l - ix1 - ix2;
    int sel = (rand() % remain);
//    if (((int)(ix1+ix2+1)*ratio)>=ix1+1) {
    if (sel<(k-ix1)) {
      train.scores[ix1].label=-1;
      train.scores[ix1].featVec=featVec;
      ++ix1;
    } else {
      test.scores[ix2].label=-1;
      test.scores[ix2].featVec=featVec;
      ++ix2;    
    }
  }
  assert(ix1==k);
  assert(ix2==l);
  while((featVec=normIter.getNext())!=NULL) {
    train.scores[ix1].label=1;
    train.scores[ix1].featVec=featVec;
    ++ix1;
    test.scores[ix2].label=1;
    test.scores[ix2].featVec=featVec;
    ++ix2;
  }
  assert(ix1==n+k);
  assert(ix2==n+l);
  train.pos=n;
  test.pos=n;
  train.neg=k;
  test.neg=l;
  train.factor = n/(double)k;
  test.factor = n/(double)l;
}

void Scores::fillFeatures(SetHandler& norm,SetHandler& shuff) {
  int n = norm.getSize()+shuff.getSize();
  ScoreHolder s;
  scores.resize(n,s);
  qVals.resize(n,-1e200); 
  int ix=0;
  const double * featVec;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  while((featVec=normIter.getNext())!=NULL) {
    scores[ix].label=1;
    scores[ix].featVec=featVec;
    ++ix;
  }
  pos=ix;
  while((featVec=shuffIter.getNext())!=NULL) {
    scores[ix].label=-1;
    scores[ix].featVec=featVec;
    ++ix;
  }
  neg=ix-pos;
  factor=pos/(double)neg;
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
    train[i].qVals.resize(train[i].scores.size(),-1e200); 
  	train[i].pos=0;train[i].neg=0;
  	for(it=train[i].begin();it!=train[i].end();it++) {
      if (it->label==1) train[i].pos++;
      else train[i].neg++;
    }
    train[i].factor=train[i].pos/(double)train[i].neg;
    test[i].qVals.resize(test[i].scores.size(),-1e200); 
    test[i].pos=0;test[i].neg=0;
  	for(it=test[i].begin();it!=test[i].end();it++) {
      if (it->label==1) test[i].pos++;
      else test[i].neg++;
  	}
    test[i].factor=test[i].pos/(double)test[i].neg;
  }
}

int Scores::calcScores(double *w,double fdr) {
  register unsigned int ix=DataSet::getNumFeatures()+1;
  if (!w_vec)
    w_vec = new double[ix];
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
  ix=0;
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
    qVals[ix++]=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (;--ix;) {
    if (qVals[ix-1]>qVals[ix]) qVals[ix-1]=qVals[ix];  
  }
  return posNow;
}

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
      if (fdr<qVals[ix1]) {
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

int Scores::getInitDirection(const double fdr, double * direction, bool findDirection) {
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


vector<double>& Scores::calcPep() {

  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));
  
  // Logistic regression on the data
  LogisticRegression lr;
  PosteriorEstimator::estimate(combined,lr,pi0);

  peps.clear();                                                                                                                  
  double old_pred = -1;
  for(size_t ix=0;ix<combined.size();++ix) {
    if (not (combined[ix].second))
      continue;
    double pred = max(lr.predict(combined[ix].first),old_pred);
    old_pred = pred;
    peps.push_back(pred);                                                                                                          }
  return peps;
}
