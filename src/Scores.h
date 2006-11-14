#ifndef SCORES_H_
#define SCORES_H_
#include <vector>
using namespace std;

class SetHandler;

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
    const vector<ScoreHolder>::const_iterator begin() const {return scores.begin();}
    const vector<ScoreHolder>::const_iterator end() const {return scores.end();}    
	int calcScores(double *w, SetHandler &set, double fdr=0.0);
//	double getPositiveTrainingIxs(const double fdr,vector<int>& set ,vector<int>& ixs);
    void getScoreAndQ(int setPos,vector<double> & s,vector<double> & fdr);
    double getQ(const double score);
	void printRoc(string & fn);	
protected:
    double *w_vec;
    const static int shortCutSize = 100;
    vector<ScoreHolder> scores;
    vector<double> qVals;
    vector<unsigned int> shortCut;
    double shortStep;
};

#endif /*SCORES_H_*/
