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
#ifndef CROSSVALIDATION_H_
#define CROSSVALIDATION_H_

#include <iostream>
#include <iomanip>
#include <vector>

#include "Globals.h"
#include "FeatureNames.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "DataSet.h"
#include "FeatureMemoryPool.h"
#include "ssl.h"

struct cpCnTriple {
  double cpos;
  double cfrac;
  unsigned int set;
  int nestedSet;
  vector<double> ww;
  int tp;
};

class CrossValidation {
  
 public:
  CrossValidation(bool quickValidation, bool reportPerformanceEachIteration, 
    double testFdr, double selectionFdr, double initialSelectionFdr, 
    double selectedCpos, double selectedCneg, int niter, bool usePi0, 
		  int nestedXvalBins, bool trainBestPositive, unsigned int numThreads);
  ~CrossValidation();
  
  int preIterationSetup(Scores & fullset, SanityCheck * pCheck, 
                        Normalizer* pNorm, FeatureMemoryPool& featurePool);
  
  void train(Normalizer* pNorm);
  
  void postIterationProcessing(Scores & fullset, SanityCheck * pCheck);
  
  void printAllWeights(ostream & weightStream, Normalizer* pNorm);
  
  void printDOC();
  void getAvgWeights(std::vector<double>& weights, Normalizer* pNorm);
  
  void inline setSelectedCpos(double cpos) { selectedCpos_ = cpos; }
  double inline getSelectedCpos() { return selectedCpos_; }
  void inline setSelectedCneg(double cneg) { selectedCneg_ = cneg; }
  double inline getSelectedCneg() { return selectedCneg_; }
  
  void inline setSelectionFdr(double fdr) { selectionFdr_ = fdr; }
  double inline getSelectionFdr() { return selectionFdr_; }
  void inline setTestFdr(double fdr) { testFdr_ = fdr; }
  
  void inline setNiter(unsigned int n) { niter_ = n; }
  unsigned int inline getNiter() { return niter_; }
  void inline setQuickValidation(bool on) { quickValidation_ = on; }
  void inline setReportPerformanceEachIteration(bool on) { 
    reportPerformanceEachIteration_ = on;
  }
  
 protected:
  std::vector<AlgIn*> svmInputs_;
  std::vector< std::vector<double> > w_; // svm weights for each fold
  std::vector<cpCnTriple> classWeightsPerFold_; // svm weights for each fold
  
  bool quickValidation_;
  bool usePi0_;
  bool reportPerformanceEachIteration_;

  unsigned int numThreads_;
  
  double testFdr_; // fdr used for cross validation performance measuring
  double selectionFdr_; // fdr used for determining positive training set
  double initialSelectionFdr_; // fdr used for determining positive training set in first iteration
  double selectedCpos_; // soft margin parameter for positive training set
  double selectedCneg_; // soft margin parameter for negative training set
  
  unsigned int niter_;
  unsigned int nestedXvalBins_;
  
  bool trainBestPositive_;
  
  const static double requiredIncreaseOver2Iterations_;
  
  const static unsigned int numFolds_;
  const static unsigned int numAlgInObjects_;
  std::vector<Scores> trainScores_, testScores_;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  void processSingleCpCnPair(cpCnTriple& cpCnFold,
			    options * pOptions, AlgIn* svmInput,
			    std::vector< std::vector<Scores> >& nestedTestScoresVec);
  void trainCpCnPair(cpCnTriple& cpCnFold,
		     options * pOptions, AlgIn* svmInput);
  void validateCpCnPairs(std::vector< std::vector<Scores> >& nestedTestScoresVec);
  int mergeCpCnPairs(vector_double* pWeights, double selectionFdr,
		     options * pOptions, std::vector< std::vector<Scores> >& nestedTestScoresVec);
  int processSingleFold(unsigned int set, double selectionFdr,
			double& best_cpos, double& best_cfrac, 
			vector_double* pWeights, options* pOptions,
			std::vector<AlgIn*> svmInputs,
			std::vector< std::vector<Scores> >& nestedTestScoresVec);
  int doStep(bool updateDOC, Normalizer* pNorm, double selectionFdr);
  
  void printSetWeights(ostream & weightStream, unsigned int set);
  void printRawSetWeights(ostream & weightStream, unsigned int set, 
                       Normalizer* pNorm);
  
  void printAllWeightsColumns(ostream & weightStream);
  void printAllRawWeightsColumns(ostream & weightStream, Normalizer* pNorm);
  static void printAllWeightsColumns(std::vector< std::vector<double> > weightMatrix, 
                              ostream & weightStream);
};

#endif /*CROSSVALIDATION_H_*/
