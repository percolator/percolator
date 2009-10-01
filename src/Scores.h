/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

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
    void fill(string & fn);
    inline unsigned int size() {return (pos+neg);}
    inline unsigned int posSize() {return (pos);}
    inline unsigned int posNowSize() {return (posNow);}
    inline unsigned int negSize() {return (neg);}
    inline static bool isOutXmlDecoys() {return outxmlDecoys;}
    inline static void setOutXmlDecoys(bool decoys_out) {outxmlDecoys=decoys_out;}
    double pi0;
    double factor;
protected:
    vector<double> w_vec;
    int neg,pos,posNow;
    double q1,q3;
    vector<ScoreHolder> scores;
    map<const double *,ScoreHolder *> scoreMap;
    DescriptionOfCorrect doc;
    static bool outxmlDecoys;
};

#endif /*SCORES_H_*/
