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
#include <cfloat>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>

#include "PSMDescription.h"
#include "FeatureNames.h"
#include "PseudoRandom.h"
#include "Normalizer.h"
#include "FeatureMemoryPool.h"

#include <boost/unordered/unordered_map.hpp>

class Scores;

/*
* ScoreHolder is a class that provides a way to assign score value to a
* PSMDescription and have a way to compare PSMs based on the assigned
* score value and output them into the stream.
*
* Here are some useful abbreviations:
* PSM - Peptide Spectrum Match
*
*/
class ScoreHolder {
 public:
  double score, q, pep, p;
  PSMDescription* pPSM;
  int label;
  
  ScoreHolder() : score(0.0), q(0.0), pep(0.0), p(0.0), label(0), pPSM(NULL) {}
  ScoreHolder(const double s, const int l, PSMDescription* psm = NULL) :
    score(s), q(0.0), pep(0.0), p(0.0), label(l), pPSM(psm) {}
  virtual ~ScoreHolder() {}
  
  std::pair<double, bool> toPair() const { 
    return pair<double, bool> (score, label > 0); 
  }
  
  inline bool isTarget() const { return label != -1; }
  inline bool isDecoy() const { return label == -1; }
  void printPSM(ostream& os, bool printDecoys, bool printExpMass);
  void printPepXML(ostream& os, map<char,float> &aaWeight, int index);
  void printPeptide(ostream& os, bool printDecoys, bool printExpMass, Scores& fullset);
  std::string getCharge(std::string id);
};
struct lessThanBaseName
{
    inline bool operator() (const ScoreHolder& struct1, const ScoreHolder& struct2)
    {
        std::string id1 = struct1.pPSM->getId();
        std::string id2 = struct2.pPSM->getId();
        return (id1 < id2);
    }
};
inline bool operator>(const ScoreHolder& one, const ScoreHolder& other);
inline bool operator<(const ScoreHolder& one, const ScoreHolder& other);
  
struct lexicOrderProb : public binary_function<ScoreHolder, ScoreHolder, bool> {
  static int compStrIt(std::string::iterator first1, std::string::iterator last1,
                       std::string::iterator first2, std::string::iterator last2) {
    for ( ; (first1 != last1) && (first2 != last2); first1++, first2++ ) {
      if (*first1 < *first2) return 1;
      if (*first2 < *first1) return -1;
    }
    if (first2 != last2) return 1;
    else if (first1 != last1) return -1;
    else return 0;
  }
  
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    int peptCmp = compStrIt(__x.pPSM->getFullPeptideSequence().begin() + 2, 
                            __x.pPSM->getFullPeptideSequence().end() - 2, 
                            __y.pPSM->getFullPeptideSequence().begin() + 2, 
                            __y.pPSM->getFullPeptideSequence().end() - 2);
    return ( ( peptCmp == 1 ) 
    || ( (peptCmp == 0) && (__x.label > __y.label) )
    || ( (peptCmp == 0) && (__x.label == __y.label) && (__x.score > __y.score) ) );
  }
};

/**
 * Computes a fast hash for unsigned int using the function suggested at:
 * https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
 * Note that std::hash is sometimes implemented as the identity function, defeating 
 * the purpose of random shuffling.
 */
inline unsigned int fast_uint_hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

/**
 * Orders ScoreHolders by a hash computed from the specFileNr and scan.
 * Uses the hash combination function h1 ^ (h2 << 1) suggested in https://en.cppreference.com/w/cpp/utility/hash
 */
struct OrderScanHash : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    size_t hash_x = fast_uint_hash(__x.pPSM->specFileNr) ^ (fast_uint_hash(__x.pPSM->scan) << 1);
    size_t hash_y = fast_uint_hash(__y.pPSM->specFileNr) ^ (fast_uint_hash(__y.pPSM->scan) << 1);
    return (hash_x < hash_y);
  }
};

struct OrderScanMassCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->specFileNr < __y.pPSM->specFileNr )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan < __y.pPSM->scan) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass < __y.pPSM->expMass) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.score > __y.score) ) );
  }
};

struct OrderScanMassLabelCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->specFileNr < __y.pPSM->specFileNr )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan < __y.pPSM->scan) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass < __y.pPSM->expMass) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.label > __y.label) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.label == __y.label) && (__x.score > __y.score) ) );
  }
};

struct OrderScanLabel : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->specFileNr < __y.pPSM->specFileNr )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan < __y.pPSM->scan) )
    || ( (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.label > __y.label) ) );
  }
};


struct UniqueScanMassCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass);
  }
};

struct UniqueScanMassLabelCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.label == __y.label) && (__x.pPSM->expMass == __y.pPSM->expMass);
  }
};

struct UniqueScanLabel : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->specFileNr == __y.pPSM->specFileNr) && (__x.pPSM->scan == __y.pPSM->scan) && (__x.label == __y.label);
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
* FDR - False Discovery Rate
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
  void merge(vector<Scores>& sv, double fdr, bool skipNormalizeScores, std::vector< std::vector<double> >& all_w);
  void postMergeStep();
  
  std::vector<ScoreHolder>::iterator begin() { return scores_.begin(); }
  std::vector<ScoreHolder>::iterator end() { return scores_.end(); }
  
  double calcScore(const double* features, const std::vector<double>& w) const;
  void scoreAndAddPSM(ScoreHolder& sh, const std::vector<double>& rawWeights,
                      FeatureMemoryPool& featurePool);
  int calcScores(vector<double>& w, double fdr, bool skipDecoysPlusOne = false);
  int calcQ(double fdr, bool skipDecoysPlusOne = false);
  void recalculateDescriptionOfCorrect(const double fdr);
  void calcPep();
  
  void populateWithPSMs(SetHandler& setHandler);
  
  int getInitDirection(const double initialSelectionFdr, std::vector<double>& direction);
  void createXvalSetsBySpectrum(std::vector<Scores>& train, 
      std::vector<Scores>& test, const unsigned int xval_fold,
      FeatureMemoryPool& featurePool);
  
  void generatePositiveTrainingSet(AlgIn& data, const double fdr,
      const double cpos, const bool trainBestPositive);
  void generateNegativeTrainingSet(AlgIn& data, const double cneg);
  
  void recalculateSizes();
  void normalizeScores(double fdr, std::vector<double>& weights);
  
  void weedOutRedundant();
  void weedOutRedundant(std::map<std::string, unsigned int>& peptideSpecCounts, 
                        double specCountQvalThreshold);
  void weedOutRedundantTDC();
  void weedOutRedundantMixMax();
  
  void printRetentionTime(ostream& outs, double fdr);
  unsigned getQvaluesBelowLevel(double level);
  
  
  void print(int label, std::ostream& os = std::cout);
    
  inline double getPi0() const { return pi0_; }
  inline double getTargetDecoySizeRatio() const { 
    return targetDecoySizeRatio_; 
  }
  inline unsigned int size() const { 
    return totalNumberOfTargets_ + totalNumberOfDecoys_; 
  }
  inline unsigned int posSize() const { return totalNumberOfTargets_; }
  inline unsigned int negSize() const { return totalNumberOfDecoys_; }  
  
  inline void addScoreHolder(const ScoreHolder& sh) {
    scores_.push_back(sh);
  }
  
  std::vector<PSMDescription*>& getPsms(PSMDescription* pPSM) {
    return peptidePsmMap_[pPSM];
  }
  
  void reset() { 
    scores_.clear(); 
    totalNumberOfTargets_ = 0;
    totalNumberOfDecoys_ = 0;
  }
  void setUsePi0(bool usePi0);
  inline void setOutputRT(const bool is_output_rt){
    is_output_rt_ = is_output_rt;
  }
 protected:
  bool usePi0_;
  
  double pi0_;
  double targetDecoySizeRatio_;
  unsigned int totalNumberOfDecoys_, totalNumberOfTargets_;
  
  std::vector<ScoreHolder> scores_;
  std::map<PSMDescription*, std::vector<PSMDescription*> > peptidePsmMap_;
  
  double* decoyPtr_;
  double* targetPtr_;
  
  void reorderFeatureRows(FeatureMemoryPool& featurePool, bool isTarget,
    boost::unordered_map<double*, double*>& movedAddresses, size_t& idx);
  void getScoreLabelPairs(std::vector<pair<double, bool> >& combined);
  void checkSeparationAndSetPi0();
  bool is_output_rt_ = false;
};

#endif /*SCORES_H_*/
