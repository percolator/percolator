/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.h,v 1.31 2007/12/04 01:49:01 lukall Exp $
 *******************************************************************************/
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
  ScoreHolder():score(0.0),featVec(NULL),label(0){;}
  virtual ~ScoreHolder() {;}
};

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other); 
inline bool operator<(const ScoreHolder &one, const ScoreHolder &other); 

class AlgIn;

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
    void static fillFeatures(Scores& train,Scores& thresh,Scores& test,SetHandler& norm,SetHandler& shuff,
                             const double trainRatio,const double testRatio);
    void createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold);
    void generatePositiveTrainingSet(AlgIn& data,const double fdr,const double cpos);
    void generateNegativeTrainingSet(AlgIn& data,const double cneg);
    int getInitDirection(const double fdr, double * direction, bool findDirection);
    double getQ(const double score);
    vector<double>& calcPep();
    void printRoc(string & fn); 
    void fill(string & fn);
    inline unsigned int size() {return (pos+neg);} 
    inline unsigned int posSize() {return (pos);} 
    inline unsigned int posNowSize() {return (posNow);} 
    inline unsigned int negSize() {return (neg);} 
    static double pi0;
    double factor;
protected:
    double *w_vec;
    int neg;
    int pos;
    int posNow;
    const static int shortCutSize = 100;
    const static int pepBins = 40;
    vector<ScoreHolder> scores;
    vector<double> peps;
    vector<double> qVals;
    vector<unsigned int> shortCut;
    double shortStep;
};

#endif /*SCORES_H_*/
