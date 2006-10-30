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

bool operator>(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score>other.score);}

bool operator<(const ScoreHolder &one, const ScoreHolder &other) 
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
  if (fdr>0.0) {
  	double tp=0,fp=0;
    vector<ScoreHolder>::iterator it,it2;
    for(it=scores.begin();it!=scores.end();it++) {
      if (it->label!=-1)
        tp++;
      if (it->label==-1)
        fp++;
      if (fdr<(fp/(tp+fp))) {
        double zero = it->score;
        w[DataSet::getNumFeatures()] -= zero;
        for(it2=scores.begin();it2!=scores.end();it2++) 
          it2->score -= zero;
        if (VERB>4) { cerr << "Scores lowered by " << zero << endl;}
        break;
      }
    }
    return (int) tp;
  }
  return 0;
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
