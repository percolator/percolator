/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.h,v 1.38 2008/06/09 16:51:45 lukall Exp $
 *******************************************************************************/
#ifndef SCORES_H_
#define SCORES_H_
#include <vector>
#include <map>
using namespace std;

class SetHandler;

class ScoreHolder{
public:
  double score,q,pep;
  const double * featVec;
  int label;
  ScoreHolder():score(0.0),q(-1.),pep(-1.),featVec(NULL),label(0){;}
  ScoreHolder(const double &s,const int &l, const double * fv = NULL):score(s),q(-1.),pep(-1.),featVec(fv),label(l){;}
  virtual ~ScoreHolder() {;}
  pair<double,bool> toPair() {return pair<double,bool>(score,label>0);}
};

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other); 
inline bool operator<(const ScoreHolder &one, const ScoreHolder &other); 

class AlgIn;

class Scores
{
public:
	Scores();
	~Scores();
    void merge(vector<Scores>& sv);
	double calcScore(const double * features) const;
    const vector<ScoreHolder>::const_iterator begin() const {return scores.begin();}
    const vector<ScoreHolder>::const_iterator end() const {return scores.end();}    
	int calcScores(vector<double>& w, double fdr=0.0);
    void fillFeatures(SetHandler& norm,SetHandler& shuff);
    void static fillFeatures(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio);
    void createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold);
    void generatePositiveTrainingSet(AlgIn& data,const double fdr,const double cpos);
    void generateNegativeTrainingSet(AlgIn& data,const double cneg);
    void normalizeScores();
    int getInitDirection(const double fdr, vector<double>& direction, bool findDirection);
 //   double getQ(const double score);
    ScoreHolder * getScoreHolder(const double *d);
    void calcPep();
    double estimatePi0();
    void printRoc(string & fn); 
    void fill(string & fn);
    inline unsigned int size() {return (pos+neg);} 
    inline unsigned int posSize() {return (pos);} 
    inline unsigned int posNowSize() {return (posNow);} 
    inline unsigned int negSize() {return (neg);} 
    static double pi0;
    double factor;
protected:
    vector<double> w_vec;
    int neg,pos,posNow;
    double q1,q3;
    vector<ScoreHolder> scores;
    map<const double *,ScoreHolder *> scoreMap; 
};

#endif /*SCORES_H_*/
