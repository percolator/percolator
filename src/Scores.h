/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the
 * Department of Genome Sciences at the University of Washington.
 *
 * $Id: Scores.h,v 1.45 2009/01/04 22:49:30 lukall Exp $
 *******************************************************************************/
#ifndef SCORES_H_
#define SCORES_H_
#include <vector>
#include <map>
#include <iostream>
using namespace std;
#include "DescriptionOfCorrect.h"
#include "PSMDescription.h"
class SetHandler;

class ScoreHolder{
public:
  double score; // ,q,pep;
  PSMDescription * pPSM;
//  const double * featVec;
  int label;
  ScoreHolder():score(0.0),pPSM(NULL),label(0){;}
  ScoreHolder(const double &s,const int &l, PSMDescription * psm = NULL):score(s),pPSM(psm),label(l){;}
  virtual ~ScoreHolder() {;}
  pair<double,bool> toPair() {return pair<double,bool>(score,label>0);}
};

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other);
inline bool operator<(const ScoreHolder &one, const ScoreHolder &other);
ostream & operator<<(ostream& os, const ScoreHolder& sh);

class AlgIn;

class Scores
{
public:
	Scores();
	~Scores();
    void merge(vector<Scores>& sv);
	double calcScore(const double * features) const;
    vector<ScoreHolder>::iterator begin() {return scores.begin();}
    vector<ScoreHolder>::iterator end() {return scores.end();}
	int calcScores(vector<double>& w, double fdr=0.0);
	int calcQ(double fdr=0.0);
    void fillFeatures(SetHandler& norm,SetHandler& shuff);
    void createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold);
    void recalculateDescriptionOfGood(const double fdr);
    void generatePositiveTrainingSet(AlgIn& data,const double fdr,const double cpos);
    void generateNegativeTrainingSet(AlgIn& data,const double cneg);
    void normalizeScores();
    void printRetentionTime(ostream& outs, double fdr);
    int getInitDirection(const double fdr, vector<double>& direction, bool findDirection);
    ScoreHolder* getScoreHolder(const double *d);
    DescriptionOfCorrect& getDOC() {return doc;}
    void setDOCFeatures();
    void calcPep();
    double estimatePi0();
    double getPi0() {return pi0;}
    void printRoc(string & fn);
    void fill(string & fn);
    inline unsigned int size() {return (pos+neg);}
    inline unsigned int posSize() {return (pos);}
    inline unsigned int posNowSize() {return (posNow);}
    inline unsigned int negSize() {return (neg);}
    double pi0;
    double factor;
protected:
    vector<double> w_vec;
    int neg,pos,posNow;
    double q1,q3;
    vector<ScoreHolder> scores;
    map<const double *,ScoreHolder *> scoreMap;
    DescriptionOfCorrect doc;
};

#endif /*SCORES_H_*/
