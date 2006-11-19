#ifndef SCORES_H_
#define SCORES_H_
#include <vector>
using namespace std;

class SetHandler;

class ScoreHolder{
public:
  double score;
  const double * featVec;
  int label;
  ScoreHolder() {score=0.0;featVec=NULL;label=0;}
  virtual ~ScoreHolder() {;}
};

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other); 
inline bool operator<(const ScoreHolder &one, const ScoreHolder &other); 

class Scores
{
public:
	Scores();
	~Scores();
	double calcScore(const double *features) const;
    const vector<ScoreHolder>::const_iterator begin() const {return scores.begin();}
    const vector<ScoreHolder>::const_iterator end() const {return scores.end();}    
	int calcScores(double *w, double fdr=0.0);
    void fillFeatures(SetHandler& norm,SetHandler& shuff);
    void static fillFeatures(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio);

//	double getPositiveTrainingIxs(const double fdr,vector<int>& set ,vector<int>& ixs);
    double getQ(const double score);
    void printRoc(string & fn); 
    void fill(string & fn); 
protected:
    double factor;
    double *w_vec;
    const static int shortCutSize = 100;
    vector<ScoreHolder> scores;
    vector<double> qVals;
    vector<unsigned int> shortCut;
    double shortStep;
};

#endif /*SCORES_H_*/
