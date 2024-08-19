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

enum class LabelType {
  DECOY = -1,
  UNDEFINED = 0,
  TARGET = 1,
  PSEUDO_TARGET = 2 // used for RESET algorithm
};

inline std::ostream& operator<<(std::ostream& os, LabelType label) {
  os << static_cast<int>(label);
  return os;
}

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
  LabelType label;
  
  ScoreHolder() : score(0.0), q(0.0), pep(0.0), p(0.0), label(LabelType::UNDEFINED), pPSM(NULL) {}
  ScoreHolder(const double s, const LabelType l, PSMDescription* psm = NULL) :
    score(s), q(0.0), pep(0.0), p(0.0), label(l), pPSM(psm) {}
  virtual ~ScoreHolder() {}
  
  std::pair<double, bool> toPair() const { 
    return pair<double, bool> (score, isTarget()); 
  }
  
  inline bool isTarget() const { return label == LabelType::TARGET || label == LabelType::PSEUDO_TARGET; }
  inline bool isDecoy() const { return label == LabelType::DECOY; }
  inline PSMDescription* getPSM() const { return pPSM; }
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
  
struct lexicOrderProb {
    static int compStrIt(const std::string& sv1, const std::string& sv2) {
        auto it1 = sv1.begin(), end1 = sv1.end();
        auto it2 = sv2.begin(), end2 = sv2.end();

        for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
            if (*it1 < *it2) return 1;
            if (*it2 < *it1) return -1;
        }
        if (it2 != end2) return 1;
        if (it1 != end1) return -1;
        return 0;
    }
    
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        const std::string& fullSeqX = x.pPSM->getFullPeptideSequence();
        const std::string& fullSeqY = y.pPSM->getFullPeptideSequence();

        std::string peptSeqX = fullSeqX.substr(2, fullSeqX.size() - 4);
        std::string peptSeqY = fullSeqY.substr(2, fullSeqY.size() - 4);

        int peptCmp = compStrIt(peptSeqX, peptSeqY);
        return ( ( peptCmp == 1 ) 
                 || ( (peptCmp == 0) && (x.label > y.label) )
                 || ( (peptCmp == 0) && (x.label == y.label) && (x.score > y.score) ) );
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
struct OrderScanHash {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        size_t hash_x = fast_uint_hash(x.pPSM->specFileNr) ^ (fast_uint_hash(x.pPSM->scan) << 1);
        size_t hash_y = fast_uint_hash(y.pPSM->specFileNr) ^ (fast_uint_hash(y.pPSM->scan) << 1);
        return hash_x < hash_y;
    }
};

struct OrderScanMassCharge {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        if (x.pPSM->specFileNr != y.pPSM->specFileNr)
            return x.pPSM->specFileNr < y.pPSM->specFileNr;
        if (x.pPSM->scan != y.pPSM->scan)
            return x.pPSM->scan < y.pPSM->scan;
        if (x.pPSM->expMass != y.pPSM->expMass)
            return x.pPSM->expMass < y.pPSM->expMass;
        return x.score > y.score;
    }
};

struct OrderScanMassLabelCharge {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        if (x.pPSM->specFileNr != y.pPSM->specFileNr)
            return x.pPSM->specFileNr < y.pPSM->specFileNr;
        if (x.pPSM->scan != y.pPSM->scan)
            return x.pPSM->scan < y.pPSM->scan;
        if (x.pPSM->expMass != y.pPSM->expMass)
            return x.pPSM->expMass < y.pPSM->expMass;
        if (x.label != y.label)
            return x.label > y.label;
        return x.score > y.score;
    }
};

struct OrderScanLabel {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        if (x.pPSM->specFileNr != y.pPSM->specFileNr)
            return x.pPSM->specFileNr < y.pPSM->specFileNr;
        if (x.pPSM->scan != y.pPSM->scan)
            return x.pPSM->scan < y.pPSM->scan;
        return x.label > y.label;
    }
};

struct UniqueScanMassCharge {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        return (x.pPSM->specFileNr == y.pPSM->specFileNr) && (x.pPSM->scan == y.pPSM->scan) && (x.pPSM->expMass == y.pPSM->expMass);
    }
};

struct UniqueScanMassLabelCharge {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        return (x.pPSM->specFileNr == y.pPSM->specFileNr) && (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label) && (x.pPSM->expMass == y.pPSM->expMass);
    }
};

struct UniqueScanLabel {
    bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
        return (x.pPSM->specFileNr == y.pPSM->specFileNr) && (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label);
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
    targetDecoySizeRatio_(1.0), nullTargetWinProb_(1.0), totalNumberOfDecoys_(0),
    totalNumberOfTargets_(0), decoyPtr_(NULL), targetPtr_(NULL) {}
  ~Scores() {}
  void merge(vector<Scores>& sv, double fdr, bool skipNormalizeScores, std::vector< std::vector<double> >& all_w);
  void postMergeStep();
  
  std::vector<ScoreHolder>::iterator begin() { return scores_.begin(); }
  std::vector<ScoreHolder>::iterator end() { return scores_.end(); }
  
  std::vector<ScoreHolder>::const_iterator begin() const { return scores_.begin(); }
  std::vector<ScoreHolder>::const_iterator end() const { return scores_.end(); }
  
  double calcScore(const double* features, const std::vector<double>& w) const;
  void scoreAndAddPSM(ScoreHolder& sh, const std::vector<double>& rawWeights,
                      FeatureMemoryPool& featurePool);
  int calcScores(vector<double>& w, double fdr, bool skipDecoysPlusOne = false);
  int onlyCalcScores(vector<double>& w);
  int calcQ(double fdr, bool skipDecoysPlusOne = false);
  void recalculateDescriptionOfCorrect(const double fdr);
  void calcPep(const bool spline = false, const bool interpol = false, const bool from_q = false);
  int calcBalancedFDR(double treshhold);
  
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
  
  
  void print(LabelType label, std::ostream& os = std::cout);
    
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
  
  inline void setNullTargetWinProb(const double nullTargetWinProb) {
    nullTargetWinProb_ = nullTargetWinProb;
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
  double targetDecoySizeRatio_, nullTargetWinProb_;
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
