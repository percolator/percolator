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

double Scores::calcScore(const double *feat) {
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

int Scores::calcScores(double *w,SetHandler & set, double fdr) {
  ScoreHolder s;
  scores.resize(set.getSize(),s);
  w_vec=w;
  int setPos=0;
  int ixPos=-1;
  const double * features;
  unsigned int ix=0;
  while((features=set.getNext(setPos,ixPos))!=NULL) {
  	scores[ix].score = calcScore(features);
  	scores[ix].label = set.getLabel(setPos);
  	scores[ix].index = ixPos;
  	scores[ix++].set = setPos;
  }
  assert(scores.size()==ix);
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
    for (ix=0;ix < 10;ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
    cerr << "10 worst scores and labels" << endl;
    for (ix=scores.size()-10;ix < scores.size();ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
  }
  qVals.resize(set.getSize()); 
  double tp=0,fp=0,q;
  int tp_at_treshold = 0;
  ix=0;
  vector<ScoreHolder>::iterator it,it2;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      tp++;
    if (it->label==-1)
      fp++;
    q=fp/(tp+fp);
    qVals[ix++]=q;
    if (fdr>0.0 && fdr<q) {
      tp_at_treshold = (int) tp;
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

void Scores::getScoreAndQ(int setPos,vector<double> & s, vector<double> & q) {
  int tp=0,fp=0;
  vector<ScoreHolder>::iterator it;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label==-1) {fp++;} else {tp++;}
    if (it->set == setPos) {
      s[it->index]=it->score;
      double fdr = fp/(double(tp)+fp);
      q[it->index]=fdr;
    }
  }
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
