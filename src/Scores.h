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

#include "LayerArithmetic.hpp"
#include "DescriptionOfCorrect.h"
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
  //std::vector<std::string> psms_list;
  
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
  void printPeptide(ostream& os, bool printDecoys, bool printExpMass, Scores& fullset);
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

struct OrderScanMassCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->scan < __y.pPSM->scan ) 
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass < __y.pPSM->expMass) )
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.score > __y.score) ) );
  }
};

struct OrderScanMassLabelCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->scan < __y.pPSM->scan ) 
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass < __y.pPSM->expMass) )
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.label > __y.label) )
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass) 
       && (__x.label == __y.label) && (__x.score > __y.score) ) );
  }
};

struct OrderScanLabel : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return ( (__x.pPSM->scan < __y.pPSM->scan ) 
    || ( (__x.pPSM->scan == __y.pPSM->scan) && (__x.label > __y.label) ) );
  }
};


struct UniqueScanMassCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->scan == __y.pPSM->scan) && (__x.pPSM->expMass == __y.pPSM->expMass);
  }
};

struct UniqueScanMassLabelCharge : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->scan == __y.pPSM->scan) && (__x.label == __y.label) && (__x.pPSM->expMass == __y.pPSM->expMass);
  }
};

struct UniqueScanLabel : public binary_function<ScoreHolder, ScoreHolder, bool> {
  bool operator()(const ScoreHolder& __x, const ScoreHolder& __y) const {
    return (__x.pPSM->scan == __y.pPSM->scan) && (__x.label == __y.label);
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
* Pi0 - prior probability of null hypothesis
* TDC - Target Decoy Competition
*
*/


class Scores {
 public:
  Scores(bool usePi0) : usePi0_(usePi0), pi0_(1.0),
    targetDecoySizeRatio_(1.0), totalNumberOfDecoys_(0),total_number_of_decoys(0),
    totalNumberOfTargets_(0), decoyPtr_(NULL), targetPtr_(NULL) , highest_fdr_calculated(0.0), fdr_has_been_calculated(false), largest_index_lt_fdr(0) {}
  ~Scores() {}
  void merge(vector<Scores>& sv, double fdr, bool skipNormalizeScores);
  void postMergeStep();
  
  std::vector<ScoreHolder>::iterator begin() { return scores_.begin(); }
  std::vector<ScoreHolder>::iterator end() { return scores_.end(); }
  
  void set_total_number_of_decoys();
  double* get_vector_of_just_decoy_scores();
  double calcScore(const double* features, const std::vector<double>& w) const;
  void scoreAndAddPSM(ScoreHolder& sh, const std::vector<double>& rawWeights,
                      FeatureMemoryPool& featurePool);
  int calcScores(vector<double>& w, double fdr, bool skipDecoysPlusOne = false);
  double get_fdr(unsigned tps, unsigned fps);
  int calcScoresLOH(const double fdr_threshold, const pair<double, bool> *const orig_combined_begin, pair<double, bool> *combined_begin, pair<double, bool> *combined_end, int num_tps_seen_so_far, int num_fps_seen_so_far, LayerArithmetic* la);
  int calcScoresSorted(vector<double>& w, double fdr, bool skipDecoysPlusOne = false, bool print_scores=false);
  int calcQ(double fdr, bool skipDecoysPlusOne = false, bool print_scores=false);
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
  void normalizeScores(double fdr);
  
  void weedOutRedundant();
  void weedOutRedundant(std::map<std::string, unsigned int>& peptideSpecCounts, 
                        double specCountQvalThreshold);
  void weedOutRedundantTDC();
  void weedOutRedundantMixMax();
  
  void printRetentionTime(ostream& outs, double fdr);
  unsigned getQvaluesBelowLevel(double level);
  
  void setDOCFeatures(Normalizer* pNorm);
  
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
  
 protected:
  bool usePi0_;

  double highest_fdr_calculated;
  bool fdr_has_been_calculated;
  int largest_index_lt_fdr;
  unsigned long total_number_of_decoys;

  double pi0_;
  double targetDecoySizeRatio_;
  int totalNumberOfDecoys_, totalNumberOfTargets_;
  
  std::vector<ScoreHolder> scores_;
  std::map<PSMDescription*, std::vector<PSMDescription*> > peptidePsmMap_;
  DescriptionOfCorrect doc_;
  
  double* decoyPtr_;
  double* targetPtr_;
  
  void reorderFeatureRows(FeatureMemoryPool& featurePool, bool isTarget,
    boost::unordered_map<double*, double*>& movedAddresses, size_t& idx);
  void getScoreLabelPairs(std::vector<pair<double, bool> >& combined);
  void checkSeparationAndSetPi0();
};

#endif /*SCORES_H_*/
