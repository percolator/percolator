/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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

#ifndef WIN32
  #include <stdint.h>
#endif

#include <vector>
#include <map>
#include <iostream>

using namespace std;
#include "DescriptionOfCorrect.h"
#include "PSMDescription.h"

/*
* ScoreHolder is a class that provides a way to assign score value to a
* PSMDescription and have a way to compare PSMs based on the assigned
* score value and output them into the stream.
*
* Here are some usefull abbreviations:
* PSM - Peptide Spectrum Match
*
*/
class ScoreHolder {
  public:
    double score, q, pep, p;
    int label;
    PSMDescription* pPSM;
    std::vector<std::string> psms_list;
    
    ScoreHolder() : score(0.0), q(0.0), pep(0.0), p(0.0), label(0), pPSM(NULL) {}
    ScoreHolder(const double& s, const int& l, PSMDescription* psm = NULL) :
      score(s), q(0.0), pep(0.0), p(0.0), label(l), pPSM(psm) {}
    virtual ~ScoreHolder() {}
    
    pair<double, bool> toPair() const { return pair<double, bool> (score, label > 0); }
    
    inline bool isTarget() const { return label != -1; }
    inline bool isDecoy() const { return label == -1; }
};
ostream& operator<<(ostream& os, const ScoreHolder& sh);

inline bool operator>(const ScoreHolder& one, const ScoreHolder& other);
inline bool operator<(const ScoreHolder& one, const ScoreHolder& other);
  
struct lexicOrderProb : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->getPeptideSequence() < __y.pPSM->getPeptideSequence() ) 
    || ( (__x.pPSM->getPeptideSequence() == __y.pPSM->getPeptideSequence()) && (__x.label > __y.label) )
    || ( (__x.pPSM->getPeptideSequence() == __y.pPSM->getPeptideSequence()) && (__x.label == __y.label)
      && (__x.score > __y.score) ) );
  }
};

struct OrderProb : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->getPeptideSequence() < __y.pPSM->getPeptideSequence() ) 
    || ( (__x.pPSM->getPeptideSequence() == __y.pPSM->getPeptideSequence()) && (__x.score > __y.score) ) );
  }
};


struct OrderScanMassCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool
  operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->scan < __y.pPSM->scan ) 
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass < __y.pPSM->expMass) )
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.score > __y.score) ) );
  }
};

inline string getRidOfUnprintablesAndUnicode(string inpString) {
  string outputs = "";
  for (unsigned int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    //NOTE signed char ranges -128 to 127
    if (((int)ch) >= 32) {
      outputs += ch;
    }
  }
  return outputs;
}

/**
 * ScoreHolder for unique peptides. Only differs in the way it is printed to 
 * a stream, therefore it is safe to create a vector slicing the object to
 * a normal ScoreHolder.
 */
class ScoreHolderPeptide : public ScoreHolder {
  public:
    ScoreHolderPeptide() : ScoreHolder() {}
    ScoreHolderPeptide(ScoreHolder& sh) : ScoreHolder(sh) {}
    ScoreHolderPeptide(const double& s, const int& l, PSMDescription* psm = NULL) :
      ScoreHolder(s, l, psm) {}
    virtual ~ScoreHolderPeptide() {}
};

// overloading output operator for class ScoreHolderPeptide
ostream& operator<<(ostream& os, const ScoreHolderPeptide& sh);

class SetHandler;
class AlgIn;

/*
* Scores is a container of ScoreHolders that allows you to do a sorted merge
* of vectors of ScoreHolder.
*
* Here are some usefull abbreviations:
* DOC - Description Of Correct
* FDR - False Discovery Rate
* LCG - Linear Congruential Generator
* Pi0 - prior probability of null hypothesis
* TDC - Target Decoy Competition
*
*/
class Scores {
  
  public:
    Scores();
    ~Scores();
    void merge(vector<Scores>& sv, double fdr=0.01, bool computePi0 = true);
    
    vector<ScoreHolder>::iterator begin() { return scores_.begin(); }
    vector<ScoreHolder>::iterator end() { return scores_.end(); }
    
    double calcScore(const double* features) const;
    int calcScores(vector<double>& w, double fdr = 0.01);
    int calcQ(double fdr = 0.01);
    void recalculateDescriptionOfCorrect(const double fdr);
    void calcPep();
    double estimatePi0();
    
    void fillFeatures(SetHandler& setHandler);
    
    int getInitDirection(const double fdr, vector<double>& direction);
    void createXvalSetsBySpectrum(std::vector<Scores>& train, 
        std::vector<Scores>& test, const unsigned int xval_fold);
    
    void generatePositiveTrainingSet(AlgIn& data, const double fdr,
        const double cpos);
    void generateNegativeTrainingSet(AlgIn& data, const double cneg);
    
    void recalculateSizes();
    void normalizeScores(double fdr=0.01);
    
    void weedOutRedundant(bool computePi0 = true);
    void weedOutRedundantTDC(bool computePi0 = true);
    
    void printRetentionTime(ostream& outs, double fdr);
    unsigned getQvaluesBelowLevel(double level);
    
    void setDOCFeatures();
    
    void fill(string& fn);
    
    ScoreHolder* getScoreHolder(const double* d);
    
    DescriptionOfCorrect& getDOC() { return doc_; }
    
    inline double getPi0() const { return pi0_; }
    inline double getTargetDecoySizeRatio() const { 
      return targetDecoySizeRatio_; 
    }
    inline unsigned int size() const { 
      return totalNumberOfTargets_ + totalNumberOfDecoys_; 
    }
    inline unsigned int posSize() const { return totalNumberOfTargets_; }
    inline unsigned int negSize() const { return totalNumberOfDecoys_; }
    
    inline static void setPrintDecoysInXml(bool decoysOut) { 
      printDecoysInXml_ = decoysOut; 
    }
    inline static bool getPrintDecoysInXml() { return printDecoysInXml_; }
    inline static void setShowExpMass(bool expmass) { showExpMass_ = expmass; }
    inline static bool getShowExpMass() { return showExpMass_; }
    
    inline void resetScoreMap() { scoreMap_.clear(); }
    
    inline static void setSeed(unsigned long s) { seed_ = s; }
    unsigned long lcg_rand();
    
    inline void addScoreHolder(const ScoreHolder& sh) {
      scores_.push_back(sh);
    }
    
  protected:
    
    vector<ScoreHolder> scores_;
    vector<double> svmWeights_;
    int totalNumberOfDecoys_, totalNumberOfTargets_;
    std::map<const double*, ScoreHolder*> scoreMap_;
    DescriptionOfCorrect doc_;
    static bool printDecoysInXml_;
    static unsigned long seed_;
    static bool showExpMass_;
    double pi0_;
    double targetDecoySizeRatio_;
};

#endif /*SCORES_H_*/
