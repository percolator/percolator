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

#ifdef WIN32
#ifndef uint32_t
#define uint32_t unsigned long
#endif
#ifndef uint64_t
#define uint64_t unsigned long long
#endif
#else

#include <stdint.h>
#endif
#include <vector>
#include <map>
#include <iostream>
using namespace std;
#include "DescriptionOfCorrect.h"
#include "PSMDescription.h"
class SetHandler;

class ScoreHolder {
  public:
    double score; // ,q,pep;
    PSMDescription* pPSM;
    //  const double * featVec;
    int label;
    ScoreHolder() :
      score(0.0), pPSM(NULL), label(0) {
      ;
    }
    ScoreHolder(const double& s, const int& l, PSMDescription* psm = NULL) :
      score(s), pPSM(psm), label(l) {
      ;
    }
    virtual ~ScoreHolder() {
      ;
    }
    pair<double, bool> toPair() {
      return pair<double, bool> (score, label > 0);
    }
    bool isTarget() {
      return label != -1;
    }
    bool isDecoy() {
      return label == -1;
    }
};

inline bool operator>(const ScoreHolder& one, const ScoreHolder& other);
inline bool operator<(const ScoreHolder& one, const ScoreHolder& other);
ostream& operator<<(ostream& os, const ScoreHolder& sh);

class AlgIn;

class Scores {
  public:
    Scores();
    ~Scores();
    void merge(vector<Scores>& sv, double fdr=0.01, bool reportUniquePeptides = false);
    double calcScore(const double* features) const;
    vector<ScoreHolder>::iterator begin() {
      return scores.begin();
    }
    vector<ScoreHolder>::iterator end() {
      return scores.end();
    }
    int calcScores(vector<double>& w, double fdr = 0.01);
    int calcQ(double fdr = 0.01);
    void fillFeatures(SetHandler& norm, SetHandler& shuff);
    void createXvalSets(vector<Scores>& train, vector<Scores>& test,
                        const unsigned int xval_fold);
    void recalculateDescriptionOfGood(const double fdr);
    void generatePositiveTrainingSet(AlgIn& data, const double fdr,
                                     const double cpos);
    void generateNegativeTrainingSet(AlgIn& data, const double cneg);
    void normalizeScores(double fdr=0.01);
    void weedOutRedundant();
    void printRetentionTime(ostream& outs, double fdr);
    int getInitDirection(const double fdr, vector<double>& direction,
                         bool findDirection);
    ScoreHolder* getScoreHolder(const double* d);
    DescriptionOfCorrect& getDOC() {
      return doc;
    }
    void setDOCFeatures();
    void calcPep();
    double estimatePi0();
    double getPi0() {
      return pi0;
    }
    void fill(string& fn);
    inline unsigned int size() {
      return (pos + neg);
    }
    inline unsigned int posSize() {
      return (pos);
    }
    inline unsigned int posNowSize() {
      return (posNow);
    }
    inline unsigned int negSize() {
      return (neg);
    }
    inline static bool isOutXmlDecoys() {
      return outxmlDecoys;
    }
    inline static void setOutXmlDecoys(bool decoys_out) {
      outxmlDecoys = decoys_out;
    }
    inline static void setSeed(uint32_t s) {
      seed = s;
    }
    uint32_t lcg_rand();
    double pi0;
    double targetDecoySizeRatio;
  protected:
    vector<double> w_vec;
    int neg, pos, posNow;
    vector<ScoreHolder> scores;
    std::map<const double*, ScoreHolder*> scoreMap;
    DescriptionOfCorrect doc;
    static bool outxmlDecoys;
    static uint32_t seed;
};

#endif /*SCORES_H_*/
