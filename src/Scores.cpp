#include<iostream>
#include<fstream>
#include<algorithm>
#include <set>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score>other.score);}

inline bool operator<(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score<other.score);}

Scores::Scores()
{
}

Scores::~Scores()
{
}

void Scores::printRoc(string & fn){
 ofstream rocStream(fn.data(),ios::out);
 vector<ScoreHolder>::iterator it;
 for(it=scores.begin();it!=scores.end();it++) {
   rocStream << (it->label==-1?-1:1) << endl;
 }
 rocStream.close();
}	

double Scores::calcScore(const double *feat) const{
  double score = 0.0;
  register int ix=0;
  for(;ix<DataSet::getNumFeatures();ix++) {
  	score += feat[ix]*w_vec[ix];
  }
  score += w_vec[ix];
//  if (w_vec[DataSet::getNumFeatures()]<0)
//    score = - score;
  return score;
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

  int pos = -1,set=0,ix1=0,ix2=0;
  const double * featVec;
  while((featVec=norm.getNext(set,pos))!=NULL) {
    if (((int)(ix1+ix2)*ratio)>ix1) {
      train.scores[ix1].label=1;
      train.scores[ix1].featVec=featVec;
      ++ix1;
    } else {
      test.scores[ix2].label=1;
      test.scores[ix2].featVec=featVec;
      ++ix2;    
    }
  }
  pos = -1,set=0;
  while((featVec=shuff.getNext(set,pos))!=NULL) {
    train.scores[ix1].label=-1;
    train.scores[ix1].featVec=featVec;
    ++ix1;
    test.scores[ix2].label=-1;
    test.scores[ix2].featVec=featVec;
    ++ix2;
  }
  assert(ix1==n+k);
  assert(ix2==n+l);
  train.factor = ix1/(n+k+l);
  test.factor = ix2/(n+k+l);
}

void Scores::fillFeatures(SetHandler& norm,SetHandler& shuff) {
  int n = norm.getSize()+shuff.getSize();
  ScoreHolder s;
  scores.resize(n,s);
  qVals.resize(n,-1e200); 
  int pos = -1,set=0,ix=0;
  const double * featVec;
  while((featVec=norm.getNext(set,pos))!=NULL) {
    scores[ix].label=1;
    scores[ix].featVec=featVec;
    ++ix;
  }
  pos = -1,set=0;
  while((featVec=shuff.getNext(set,pos))!=NULL) {
    scores[ix].label=-1;
    scores[ix].featVec=featVec;
    ++ix;
  }
  factor=1.0;
}


void Scores::createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold) {
 vector<vector<DataSet *> > minors(subsets.size()*xval_fold);
  for(unsigned int j=0;j<xval_fold;j++) {
    minors[j].resize(subsets.size());
    for(unsigned int i=0;i<subsets.size();i++) {
      minors[j][i]=new DataSet(*(subsets[i]),xval_fold,j);
    }  
  }
  for(unsigned int j=0;j<xval_fold;j++) {
    vector<DataSet *> ff(0);
    for(unsigned int i=0;i<xval_fold;i++) {
      if (i==j)
        continue;
      for(unsigned int k=0;k<subsets.size();k++) {
        ff.push_back(new DataSet(*minors[i][k]));
      }
    }
    train[j].setSet(ff);
    test[j].setSet(minors[j]);
  }  
}

int Scores::calcScores(double *w,double fdr) {
  w_vec=w;
  const double * features;
  vector<ScoreHolder>::iterator it = scores.begin(),it2;
  while(it!=scores.end()) {
    features = it->featVec;
  	it->score = calcScore(features);
    it++;
  }
  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
/*  double lastScore = 1e200;
  unsigned int p=0,n=0;
  for (ix=0;ix<scores.size();ix++) {
    if (scores[ix].score>lastScore)
      cerr << "ajajaj " << ix << endl;
    lastScore = scores[ix].score;
    if (scores[ix].label==1) {p++;}
    else if (scores[ix].label==-1) {n++;}
    else cerr << "Strange label " << scores[ix].label << endl;
  }
  cerr << p << " +1 labels and " << n << " -1 labels" << endl;
*/
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
  int tp=0,fp=0;
  int tp_at_treshold = 0;
  double scaled_fp,q;
  unsigned int ix=0;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      tp++;
    if (it->label==-1) {
      fp++;
      scaled_fp=fp/factor;
    }
    q=scaled_fp/(tp+scaled_fp);
    qVals[ix++]=q;
    if (fdr>0.0 && fdr<q) {
      tp_at_treshold = tp;
      double zero = it->score;
      w[DataSet::getNumFeatures()] -= zero;
      for(it2=scores.begin();it2!=scores.end();it2++) 
        it2->score -= zero;
      if (VERB>4) { cerr << "Scores lowered by " << zero << endl;}
      fdr = -1;
    }
  }
  return tp_at_treshold;
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

/*
double Scores::getPositiveTrainingIxs(const double fdr,vector<int>& set,vector<int>& ixs) {
  double tp=0,fp=0;
  vector<ScoreHolder>::iterator it;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label==-1) fp++; else tp++;
    if (fdr<(fp/(tp+fp)))
      return -it->score;
    if (it->label!=-1)
      ixs.push_back(it->index);
      set.push_back(it->set);      
  }
  return 0;
}
*/
