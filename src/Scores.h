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
// #ifndef uint32_t
// #define uint32_t unsigned long
// #endif
// #ifndef uint64_t
// #define uint64_t unsigned long long
// #endif
#else

#include <stdint.h>
#endif
#include <vector>
#include <map>
#include <iostream>
using namespace std;
#include "DescriptionOfCorrect.h"
#include "PSMDescription.h"

#include "percolator_out.hxx"

class SetHandler;

class ScoreHolder {
  public:
    double score; // ,q,pep;
    PSMDescription* pPSM;
    //  const double * featVec;
    int label;
    vector<string> psms_list;
    ScoreHolder() :
      score(0.0), label(0), pPSM(NULL), psms_list () {
      ;
    }
    ScoreHolder(const double& s, const int& l, PSMDescription* psm = NULL) :
      score(s), label(l), pPSM(psm), psms_list () {
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
std::auto_ptr< ::percolatorOutNs::psm> returnXml_PSM(const vector<ScoreHolder>::iterator);
ostream& operator<<(ostream& os, const ScoreHolder& sh);

struct lexicOrder : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return __x.pPSM->getPeptide() < __y.pPSM->getFullPeptide();
  }
};

struct lexicOrderProb : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->getPeptide() < __y.pPSM->getFullPeptide() ) 
    || ( (__x.pPSM->getPeptide() == __y.pPSM->getFullPeptide()) && (__x.label != __y.label) )
    || ( (__x.pPSM->getPeptide() == __y.pPSM->getFullPeptide()) && (__x.label == __y.label)
      && (__x.score > __y.score) ) );
  }
};

struct lexicEq : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    string xPept = __x.pPSM->getPeptideSequence();
    string yPept = __y.pPSM->getPeptideSequence();
    if(xPept.compare(yPept) == 0) return true;
    else return false;
  }
};

inline string getRidOfUnprintablesAndUnicode(string inpString) {
  string outputs = "";
  for (int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    if (((int)ch) >= 32 && ((int)ch) <= 128) {
      outputs += ch;
    }
  }
  return outputs;
}

/**
 * ScoreHolder for unique peptides.
 */
class ScoreHolderPeptide: public ScoreHolder {
  public:
    ScoreHolderPeptide() :
      ScoreHolder(){
      ;
    }
    ScoreHolderPeptide(ScoreHolder& sh) :
      ScoreHolder(sh){
      ;
    }
    ScoreHolderPeptide(const double& s, const int& l, PSMDescription* psm = NULL) :
      ScoreHolder(s, l, psm) {
    	;
    }
    virtual ~ScoreHolderPeptide() {
      ;
    }
};

// overloading output operator for class ScoreHolderPeptide
ostream& operator<<(ostream& os, const ScoreHolderPeptide& sh);


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
    void fillFeatures(SetHandler& norm, SetHandler& shuff, bool);
    void createXvalSets(vector<Scores>& train, vector<Scores>& test,
        const unsigned int xval_fold);
    void createXvalSetsBySpectrum(vector<Scores>& train, vector<Scores>& test,
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
    
    /** Return a list of peptides that are matching the given protein name **/
    std::vector<std::string> getPeptides(std::string proteinName);
    
    /** Return the scores whose q value is less or equal than the threshold given**/
    unsigned getQvaluesBelowLevel(double level);
    
    void fill(string& fn);
    inline unsigned int size() {
      return (totalNumberOfTargets + totalNumberOfDecoys);
    }
    inline unsigned int posSize() {
      return (totalNumberOfTargets);
    }
    inline unsigned int posNowSize() {
      return (posNow);
    }
    inline unsigned int negSize() {
      return (totalNumberOfDecoys);
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
    inline static void setShowExpMass(bool expmass) {
      showExpMass = expmass;
    }
    inline static bool getShowExpMass() {
      return showExpMass;
    }
    
    uint32_t lcg_rand();
    double pi0;
    double targetDecoySizeRatio;
    vector<ScoreHolder> scores;
  protected:
    vector<double> w_vec;
    int totalNumberOfDecoys, totalNumberOfTargets, posNow;
    std::map<const double*, ScoreHolder*> scoreMap;
    DescriptionOfCorrect doc;
    static bool outxmlDecoys;
    static uint32_t seed;
    static bool showExpMass;
};

#endif /*SCORES_H_*/
