#include<iostream>
#include<fstream>
#include<algorithm>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"

bool operator>(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score>other.score);}

bool operator<(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score<other.score);}

Scores::Scores()
{
}

Scores::~Scores()
{
/*	vector<ScoreHolder>::iterator it;
	it=scores.begin();
	while(it!=scores.end()) {
		delete &(*it);
		*it=NULL;
		it++;
    } 
*/
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
  for(int ix=0;ix<NUM_FEATURE;ix++) {
  	score += feat[ix]*w_vec[ix];
  }
//  if (w_vec[NUM_FEATURE]<0)
//    score = - score;
  return score;
}

void Scores::calcScores(double *w,IsoChargeSet & set) {
  ScoreHolder s;
  scores.resize(set.getSize(),s);
  w_vec=w;
  int setPos=0;
  int ixPos=-1;
  const double * features;
  int ix=-1;
  while((features=set.getNext(setPos,ixPos))!=NULL) {
  	scores[++ix].score = calcScore(features);
  	scores[ix].label = set.getLabel(&setPos);
  	scores[ix].index = ixPos;
  	scores[ix].set = setPos;
  }
  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
/*  for (ix=0;ix < 10;ix++) {
  	cout << scores[ix].score << " " << scores[ix].label << endl;
  }
  for (ix=scores.size()-10;ix < scores.size();ix++) {
  	cout << scores[ix].score << " " << scores[ix].label << endl;
  } */
}

double Scores::getPositiveTrainingIxs(const double fdr,vector<int>& ixs) {
  double tp=0,fp=0;
  vector<ScoreHolder>::iterator it;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      tp++;
    if (it->label==-1)
      fp++;
    if (fdr<(fp/(tp+fp)))
      return -it->score;
    if (it->label!=-1)
      ixs.push_back(it->index);
  }
  return 0;
}
