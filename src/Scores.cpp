/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.cpp,v 1.75 2008/08/27 13:57:04 lukall Exp $
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
  for(vector<Scores>::iterator a = sv.begin();a!=sv.end();a++) {
  	a->normalizeScores();
    copy(a->begin(),a->end(),back_inserter(scores));
  }
  sort(scores.begin(),scores.end(), greater<ScoreHolder>() );
}

void Scores::printRetentionTime(ostream& outs, double fdr){
  vector<ScoreHolder>::iterator it;
  for(it=scores.begin();it!=scores.end() && it->pPSM->q <= fdr; ++it) {
    if (it->label!=-1)
      outs << it->pPSM->retentionTime << "\t" << doc.estimateRT(it->pPSM->retentionFeatures) << "\t" << it->pPSM->peptide << endl;
  }
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
  register int ix=FeatureNames::getNumFeatures();
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
      scoreMap[it->pPSM->features] = &(*it);
    }
  }
  return scoreMap[d];  
}


void Scores::fillFeatures(SetHandler& norm,SetHandler& shuff) {
  scores.clear();
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  while((pPSM=normIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,1,pPSM));
  }
  while((pPSM=shuffIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,-1,pPSM));
  }
  pos = norm.getSize(); neg = shuff.getSize();
  factor=norm.getSize()/(double)shuff.getSize();
}


void Scores::createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold) {
  train.resize(xval_fold);
  test.resize(xval_fold);
  vector<size_t> remain(xval_fold);
  size_t fold = xval_fold, ix = scores.size();
  while (fold--) {
    remain[fold] = ix / (fold + 1);
    ix -= remain[fold];
  }
  
  for(unsigned int j=0;j<scores.size();j++) {
    ix = rand()%(scores.size()-j);
    fold = 0;
    while(ix>remain[fold])
      ix-= remain[fold++];
    for(unsigned int i=0;i<xval_fold;i++) {
      if(i==fold) {
        test[i].scores.push_back(scores[j]);
      } else {
        train[i].scores.push_back(scores[j]);
      }
    }
    --remain[fold];  
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

void Scores::normalizeScores() {
  // Normalize scores so that distance between 1st and 3rd quantile of the null scores are 1
  unsigned int q1index = neg/4,q3index = neg*3/4,decoys=0;
  vector<ScoreHolder>::iterator it = scores.begin();
  double q1 = it->score;
  double q3 = q1 + 1.0;
  while(it!=scores.end()) {
  	if (it->label == -1) {
  	  if(++decoys==q1index)
  	    q1=it->score;
  	  else if (decoys==q3index) {
  	    q3=it->score;
  	    break;
  	  }
  	}
    ++it;
  }
  double diff = q1-q3;
  assert(diff>0);
  for(it=scores.begin();it!=scores.end();++it) { 
    it->score -= q1;
    it->score /= diff;
  }
}

int Scores::calcScores(vector<double>& w,double fdr) {
  w_vec=w;
  const double * features;
  unsigned int ix;
  vector<ScoreHolder>::iterator it = scores.begin();
  while(it!=scores.end()) {
    features = it->pPSM->features;
  	it->score = calcScore(features);
    it++;
  }
  sort(scores.begin(),scores.end(),greater<ScoreHolder>());
  if (VERB>3) {
    cerr << "10 best scores and labels" << endl;
    for (ix=0;ix < 10;ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
    cerr << "10 worst scores and labels" << endl;
    for (ix=scores.size()-10;ix < scores.size();ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
  }
  return calcQ(fdr);
}

int Scores::calcQ(double fdr) {
  vector<ScoreHolder>::iterator it;
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
    it->pPSM->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (int ix=scores.size();--ix;) {
    if (scores[ix-1].pPSM->q > scores[ix].pPSM->q)
      scores[ix-1].pPSM->q = scores[ix].pPSM->q;  
  }
  return posNow;
}

void Scores::generateNegativeTrainingSet(AlgIn& data,const double cneg) {
  unsigned int ix1=0,ix2=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==-1) {
      data.vals[ix2]=scores[ix1].pPSM->features;
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
      if (fdr<scores[ix1].pPSM->q) {
        posNow=p;
        break;
      }
      data.vals[ix2]=scores[ix1].pPSM->features;
      data.Y[ix2]=1;
      data.C[ix2++]=cpos;
      ++p;
    }
  }
  data.positives=p;
  data.m=ix2;
}

void Scores::recalculateDescriptionOfGood(const double fdr) {
  doc.clear();
  unsigned int ix1=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==1) {
      if (fdr>scores[ix1].pPSM->q) {
        doc.registerCorrect(scores[ix1].pPSM);
      }
    }
  }
  doc.trainCorrect();
  setDOCFeatures();
}

void Scores::setDOCFeatures() {
  for(size_t ix1=0;ix1<size();++ix1) {
    doc.setFeatures(scores[ix1].pPSM);
  }
}


int Scores::getInitDirection(const double fdr, vector<double>& direction, bool findDirection) {
  int bestPositives = -1;
  int bestFeature =-1;
  bool lowBest = false;
  
  if (findDirection) { 
    for (unsigned int featNo=0;featNo<FeatureNames::getNumFeatures();featNo++) {
      vector<ScoreHolder>::iterator it = scores.begin();
      while(it!=scores.end()) {
        it->score = it->pPSM->features[featNo];
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
    for (int ix=FeatureNames::getNumFeatures();ix--;) {
      direction[ix]=0;
    }
    direction[bestFeature]=(lowBest?-1:1);
    if (VERB>1) {
      cerr << "Selected feature number " << bestFeature +1 << " as initial search direction, could separate " << 
              bestPositives << " positives in that direction" << endl;
    }
  } else {
    bestPositives = calcScores(direction,fdr);
  }
  return bestPositives;
}

double Scores::estimatePi0() {
  vector<pair<double,bool> > combined;
  vector<double> pvals;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));

  // Estimate pi0
  PosteriorEstimator::getPValues(combined,pvals);
  pi0 = PosteriorEstimator::estimatePi0(pvals);
  return pi0;

}

void Scores::calcPep() {

  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));
  vector<double> peps;                                                                                                                  
  
  // Logistic regression on the data
  PosteriorEstimator::estimatePEP(combined,pi0,peps,true);

  for(size_t ix=0; ix<scores.size(); ix++) {
    (scores[ix]).pPSM->pep = peps[ix];
  }
}
