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
 
#ifndef PROTEIN_PROB_ESTIMATOR_H_
#define PROTEIN_PROB_ESTIMATOR_H_

#include <functional>
#include <numeric>
#include <iterator>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Globals.h"
#include "ProteinFDRestimator.h"
#include "ProteinScoreHolder.h"
#include "PosteriorEstimator.h"
#include "Scores.h"
#include "PseudoRandom.h"

/** set of helper functions to sort data structures and some operations overloaded **/
struct IntCmpProb {
  bool operator()(const ProteinScoreHolder& lhs, const ProteinScoreHolder& rhs) {
    return (  (lhs.getPEP() < rhs.getPEP())
         || ( (lhs.getPEP() == rhs.getPEP()) && (lhs.getGroupId() < rhs.getGroupId()) )
         || ( (lhs.getPEP() == rhs.getPEP()) && (lhs.getGroupId() == rhs.getGroupId())
            && (lhs.getQ() < rhs.getQ()) )
         || ( (lhs.getPEP() == rhs.getPEP()) && (lhs.getGroupId() == rhs.getGroupId())
            && (lhs.getQ() == rhs.getQ()) && (lhs.getName() < rhs.getName()) )  
    );
  }
};

struct IntCmpScore {
  bool operator()(const ProteinScoreHolder& lhs, const ProteinScoreHolder& rhs) {
    return ( lhs.getScore() < rhs.getScore()
         || (lhs.getScore() == rhs.getScore() && lhs.getGroupId() < rhs.getGroupId())
         || (lhs.getScore() == rhs.getScore() && lhs.getGroupId() == rhs.getGroupId()
            && lhs.getName() < rhs.getName())
    );
  }
};
  
inline double myminfunc(double a, double b) {
  return a > b ? b : a;
}

struct RetrieveKey {
  template <typename T>
  typename T::first_type operator()(T keyValuePair) const {
    return keyValuePair.first;
  }
};
    
struct RetrieveValue {
  template <typename T>
  typename T::second_type operator()(T keyValuePair) const {
    return keyValuePair.second;
  }
};

/*
* ProteinProbEstimator is a class that computes probabilities and statistics based
* on provided proteins from the set of scored peptides from percolator.
*
* Here are some usefull abbreviations:
* Mayu - a software package for the analysis of (large) mass
*           spectrometry-based shotgun proteomics data sets.
*/
class ProteinProbEstimator { 
 public:  
  /** PROTEIN FDR ESTIMATOR PARAMETERS **/
  
  /* Default configuration (changeable by functions)
   * decoy prefix = random
   * number of bins = 10
   * target decoy ratio = 1.0
   * binning mode = equal deepth
   * correct identical sequences = true
   */
  
  const static bool correct_identical_sequences = true;
  const static bool binning_equal_deepth = true;
  const static double target_decoy_ratio;
  const static unsigned number_bins = 10;
  
  
  /** GENERAL PARAMETERS **/
  
  /** compute global protein FDR using MAYU based implementation **/
  const static bool mayufdr = false;
  /** threshold used to estimate the protein FDR(pi0) **/
  const static double psmThresholdMayu;
  /** whether to use the decoy prefix(faster) to check if a protein is decoy or not **/
  const static bool useDecoyPrefix = true;
  /** whether to count decoy proteins when estimated q values or not **/
  const static bool countDecoyQvalue_ = true;
  /** protein prior probability used to estimate the peptide prior probabilities **/
  const static double prior_protein;
  
  
  /******************************************************************************************************************/
  
  
  ProteinProbEstimator(bool trivialGrouping = true, double absenceRatio = 1.0, 
      bool outputEmpirQVal = false, std::string decoyPattern = "random");
  
  virtual ~ProteinProbEstimator();
  
  /** reads the proteins from the set of scored peptides from percolator **/
  virtual bool initialize(Scores& peptideScores);
  
  /** start the protein probabilities tool**/
  virtual void run() = 0;
  
  /** initialize the estimation of the protein probabilities **/
  virtual void computeProbabilities(const std::string& fname = "") = 0;
  
  /** print out the tab delimited list of proteins to std::cerr or the screen */
  void printOut(const std::string &proteinFN, 
     const std::string &proteinDecoyFN);
  
  /** initialize the estimation of the q values and p values **/
  void computeStatistics();
  
  /** MAYUS method for estimation of Protein FDR **/
  void computeFDR();
  
  /** write the list of proteins to the output file **/
  void writeOutputToXML(string xmlOutputFN, bool outputDecoys);

  /** Return the number of proteins whose q value is less or equal than the threshold given**/
  unsigned getQvaluesBelowLevel(double level);
  unsigned getQvaluesBelowLevelDecoy(double level);
 
  /** populate the list of proteins**/
  void setTargetandDecoysNames(Scores& peptideScores);
  
  /** add proteins read from the database **/
  void addProteinDb(bool isDecoy, std::string name, std::string sequence, double length);
  
  /** print copyright of the author**/
  virtual string printCopyright() = 0;
  
  /** printing custom parameters to xml **/
  virtual std::ostream& printParametersXML(std::ostream &os) = 0;
  
  /**some getters and setters**/
  static inline void setCalcProteinLevelProb(bool on) { calcProteinLevelProb = on; }
  static inline bool getCalcProteinLevelProb() { return calcProteinLevelProb; }
  inline bool getUsePi0() { return usePi0_; }
  double getPi0() { return pi0_; }
  double getAbsenceRatio() { return absenceRatio_; }
  double getFDR() { return fdr_; }
  
  /** return the data structure for the proteins. DO NOT REMOVE, used by Crux! **/
  const std::vector<ProteinScoreHolder>& getProteinsByRef() const { return proteins_; }
  std::vector<ProteinScoreHolder> getProteins() const { return proteins_; }
  
 protected:

  static bool calcProteinLevelProb;
  
  inline bool lastProteinInGroup(
      std::vector<ProteinScoreHolder>::const_iterator it) {
    return !trivialGrouping_ || it+1 != proteins_.end() 
                             || it->getGroupId() != (it+1)->getGroupId();
  }
  /** functions to count number of target and decoy proteins **/
  unsigned countTargets(const std::vector<std::string> &proteinList);
  unsigned countDecoys(const std::vector<std::string> &proteinList);
  bool isTarget(const std::string& proteinName);
  bool isDecoy(const std::string& proteinName);
  
   /** print a tab delimited list of proteins probabilities in a file or stdout**/
  void print(ostream& myout, bool decoy=false);
  
  /** function that extracts a list of proteins from the peptides that have a qvalue lower than psmThresholdMayu
   * this function is used to estimate the protein FDR**/
  void getTPandPFfromPeptides(double threshold, std::set<std::string> &numberTP, 
        std::set<std::string> &numberFP);
  
  /** this function generates a vector of pair protein pep and label **/
  void getCombinedList(std::vector<std::pair<double , bool> >& combined);
  
  /** compute estimated qvalues from the PEP**/
  void estimateQValues();
  
  /** compute pvalues from the scored target/decoy proteins**/
  void estimatePValues(); 
  
  /** compute empirical qvalues from the target/decoy proteins**/
  void estimateQValuesEmp();
  
  /** compute pi0 from the set of pvalues**/
  double estimatePi0(const unsigned int numBoot = 100);
  
  
  /** variables **/
  
  /** contains all the protein names for target and decoy set respectively **/
  std::set<string> truePosSet_, falsePosSet_;
  
  /** map from protein name to its sequence and its sequence length, used for Mayu method **/
  std::map<std::string,std::pair<std::string,double> > targetProteins_;
  std::map<std::string,std::pair<std::string,double> > decoyProteins_;
  
  /** vector of protein scores **/
  std::vector<ProteinScoreHolder> proteins_;
  std::map<std::string, size_t> proteinToIdxMap_;
  
  /* protein groups are either present or absent and cannot be partially present */
  bool trivialGrouping_;
  
  bool usePi0_;
  double pi0_, absenceRatio_;
  bool outputEmpirQVal_;
  double fdr_;
  unsigned int numberDecoyProteins_;
  unsigned int numberTargetProteins_;
  std::string decoyPattern_;
};

#endif /* PROTEIN_PROB_ESTIMATOR_H_ */
