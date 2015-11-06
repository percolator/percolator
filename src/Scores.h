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

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>

#include "DescriptionOfCorrect.h"
#include "PSMDescription.h"
#include "FeatureNames.h"

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
  
  std::pair<double, bool> toPair() const { 
    return pair<double, bool> (score, label > 0); 
  }
  
  inline bool isTarget() const { return label != -1; }
  inline bool isDecoy() const { return label == -1; }
  void printPSM(ostream& os, bool printDecoys, bool printExpMass);
  void printPeptide(ostream& os, bool printDecoys, bool printExpMass);
};

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
  Scores(bool usePi0) : usePi0_(usePi0), pi0_(1.0), 
    targetDecoySizeRatio_(1.0), totalNumberOfDecoys_(0),
    totalNumberOfTargets_(0), decoyPtr_(NULL), targetPtr_(NULL) {}
  ~Scores() {}
  void merge(vector<Scores>& sv, double fdr);
  void postMergeStep();
  
  std::vector<ScoreHolder>::iterator begin() { return scores_.begin(); }
  std::vector<ScoreHolder>::iterator end() { return scores_.end(); }
  
  double calcScore(const double* features, const std::vector<double>& w) const;
  int calcScores(vector<double>& w, double fdr);
  int calcQ(double fdr);
  void recalculateDescriptionOfCorrect(const double fdr);
  void calcPep();
  void estimatePi0();
  
  void fillFeatures(SetHandler& setHandler);
  
  int getInitDirection(const double fdr, std::vector<double>& direction);
  void createXvalSetsBySpectrum(std::vector<Scores>& train, 
      std::vector<Scores>& test, const unsigned int xval_fold);
  
  void generatePositiveTrainingSet(AlgIn& data, const double fdr,
      const double cpos);
  void generateNegativeTrainingSet(AlgIn& data, const double cneg);
  
  void recalculateSizes();
  void normalizeScores(double fdr);
  
  void weedOutRedundant();
  void weedOutRedundantTDC();
  
  void deleteContiguousMemoryBlock();
  
  void printRetentionTime(ostream& outs, double fdr);
  unsigned getQvaluesBelowLevel(double level);
  
  void setDOCFeatures();
  
  void print(int label, std::ostream& os = std::cout);
  
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
  
  inline static void setSeed(unsigned long s) { seed_ = s; }
  static unsigned long lcg_rand();
  
  inline void addScoreHolder(const ScoreHolder& sh) {
    scores_.push_back(sh);
  }
  
  void reset() { scores_.clear(); }
  
 protected:
  static unsigned long seed_;
  bool usePi0_;
  
  double pi0_;
  double targetDecoySizeRatio_;
  int totalNumberOfDecoys_, totalNumberOfTargets_;
  
  std::vector<ScoreHolder> scores_;
  DescriptionOfCorrect doc_;
  
  double* decoyPtr_;
  double* targetPtr_;
  
  void copyIntoContiguousMemoryBlock();
};

#endif /*SCORES_H_*/
