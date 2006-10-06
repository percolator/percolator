#include<algorithm>
#include "Scores.h"

Scores::Scores()
{
}

Scores::~Scores()
{
	vector<ScoreHolder>::iterator it;
	it=scores.begin();
	while(it!=scores.end()) {
		delete &(*it);
		it++;
    } 
}

double Scores::calcScore(const double *feat) {
  double score = 0.0;
  for(int ix=0;ix<NUM_FEATURE;ix++) {
  	score += feat[ix]*w_vec[ix];
  }
  if (w_vec[NUM_FEATURE]<0)
    score = - score;
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
  reverse(scores.begin(),scores.end());
}

void Scores::getPositiveTrainingIxs(double fdr,vector<int>& ixs) {
  double tp=0,tn=0;
  vector<ScoreHolder>::iterator it;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      tp++;
    if (it->label==-1)
      tn++;
    if (fdr<(tp/(tp+tn)))
      break;
    if (it->label!=-1)
      ixs.push_back(it->index);
  }
}
