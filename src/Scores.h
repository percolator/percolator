#ifndef SCORES_H_
#define SCORES_H_
#include <vector>
using namespace std;

class ScoreHolder{
public:
  double score;
  int index;
  int label;
  int set;
  ScoreHolder() {score=0.0;index=0;label=0;set=0;}
  ScoreHolder(const ScoreHolder& other)
    {score=other.score;index=other.index;label=other.label;set=other.set;}
  virtual ~ScoreHolder() {;}
  void operator = (const ScoreHolder& other)
    {score=other.score;index=other.index;label=other.label;set=other.set;}
};


class Scores
{
public:
	Scores();
	~Scores();
	double calcScore(const double *features);
	void calcScores(double *w, SetHandler &set);
	double getPositiveTrainingIxs(const double fdr,vector<int>& ixs);
    void Scores::getScoreAndFdr(int setPos,vector<double> & s,vector<double> & fdr);
	void printRoc(string & fn);	
protected:
    double *w_vec;
    vector<ScoreHolder> scores;
};

#endif /*SCORES_H_*/
