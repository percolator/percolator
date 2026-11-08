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

#include <iomanip>
#include <iostream>
#include <vector>

#include "DataSet.h"
#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "Globals.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "ssl.h"

struct CandidateCposCfrac {
  double cpos;
  double cfrac;
  unsigned int set;
  int nestedSet;
  vector<double> ww;
  int tp;  // validation count; 0 for the only trained pair of a fold
};

class CrossValidation {
 public:
  CrossValidation(bool quickValidation,
                  bool reportPerformanceEachIteration,
                  double testFdr,
                  double selectionFdr,
                  double initialSelectionFdr,
                  double selectedCpos,
                  double selectedCneg,
                  unsigned int niter,
                  bool usePi0,
                  unsigned int nestedXvalBins,
                  bool trainBestPositive,
                  unsigned int numThreads,
                  bool skipNormalizeScores,
                  double decoyFractionTraining,
                  unsigned int numFolds);

  int preIterationSetup(Scores& fullset,
                        SanityCheck& check,
                        const Normalizer* pNorm);

  void train(const Normalizer* pNorm);

  void postIterationProcessing(Scores& fullset, SanityCheck& check);

  void printAllWeights(ostream& weightStream, const Normalizer* pNorm);

  void getAvgWeights(std::vector<double>& weights, const Normalizer* pNorm);
  // the weights of each fold for unnormalized features, as getAvgWeights
  void getRawFoldWeights(std::vector<std::vector<double> >& weights,
                         const Normalizer* pNorm);
  // (input line, fold) of every PSM in the test sets, sorted by line, as
  // recorded by postIterationProcessing when setRecordTestFolds(true) was
  // called; the test sets themselves are emptied by the merge. Empty if
  // the sanity check rejected the trained weights. Moves the list out.
  void getTestFolds(
      std::vector<std::pair<unsigned int, unsigned int> >& lineFolds);

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
  // train as Percolator 3.09 did, see Caller's --classic-training
  void inline setClassicTraining(bool on) { classicTraining_ = on; }
  // train the SVM on the feature rows in place instead of on copies (see
  // AlgIn), which needs rows that end with the intercept 1.0, as those of a
  // FeatureMemoryPool for which canTrainInPlace(numFeatures) holds; set
  // before preIterationSetup
  void inline setSvmRowsInPlace(bool on) { svmRowsInPlace_ = on; }
  // let postIterationProcessing list the input line and fold of every test
  // PSM for getTestFolds, before the merge empties the test sets (-N)
  void inline setRecordTestFolds(bool on) { recordTestFolds_ = on; }
  // how the scores of the test bins are put on a common scale before they
  // are merged, see Caller's --merge-normalization
  void inline setMergeNormalization(MergeNormalization scheme) {
    mergeNormalization_ = scheme;
    usedMergeNormalization_ = scheme;
  }
  // the scheme postIterationProcessing merged with: the one set above, or
  // the scheme of the direction check's second chance when that kept the
  // trained SVMs (SanityCheck::secondChanceMergeScheme)
  MergeNormalization usedMergeNormalization() const {
    return usedMergeNormalization_;
  }
  // The scheme for numFolds folds when none is chosen: test up to three
  // folds, the default, whose test bins hold a third of the PSMs or more;
  // train-decoys with more folds, since the q = fdr anchor of test is noisy
  // in smaller test bins (see Caller's --num-folds).
  static MergeNormalization defaultMergeNormalization(unsigned int numFolds);

 protected:
  std::vector<AlgIn> svmInputs_;
  std::vector<std::vector<double> > weights_;  // svm weights for each fold
  std::vector<CandidateCposCfrac>
      classWeightsPerFold_;  // cpos, cneg pairs to train for each nested CV
                             // fold

  bool quickValidation_;
  bool classicTraining_ = false;
  bool svmRowsInPlace_ = false;
  bool recordTestFolds_ = false;
  std::vector<std::pair<unsigned int, unsigned int> > testFolds_;
  bool usePi0_;
  bool reportPerformanceEachIteration_;

  unsigned int numThreads_;

  double testFdr_;       // fdr used for cross validation performance measuring
  double selectionFdr_;  // fdr used for determining positive training set
  double initialSelectionFdr_;  // fdr used for determining positive training
                                // set in first iteration
  double selectedCpos_;  // soft margin parameter for positive training set
  double selectedCneg_;  // soft margin parameter for negative training set

  unsigned int niter_;
  unsigned int nestedXvalBins_;

  bool trainBestPositive_;
  bool skipNormalizeScores_;
  MergeNormalization mergeNormalization_ = MergeNormalization::TEST;
  MergeNormalization usedMergeNormalization_ = MergeNormalization::TEST;

  double decoyFractionTraining_;

  const static double requiredIncreaseOver2Iterations_;
  const static double requiredIncreaseOver2IterationsDefault_;
  const static unsigned int gridSearchIterations_;

  // after gridSearchIterations_ only each fold's selected (cpos, cfrac)
  // pair is trained
  unsigned int iteration_ = 0;
  bool gridFixed_ = false;
  std::vector<std::pair<double, double> > selectedPairs_;
  bool isSelectedPair(const CandidateCposCfrac& pair) const {
    return !gridFixed_ || (pair.cpos == selectedPairs_[pair.set].first &&
                           pair.cfrac == selectedPairs_[pair.set].second);
  }

  unsigned int numFolds_;  // number of folds for cross validation
  std::vector<Scores> trainScores_, testScores_;
  // whether trainScores_ are scored and sorted with the current weights_
  bool trainScoresCurrent_ = false;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  // buffers of the SVM solves of one thread, reused for its pairs
  struct SvmTrainingBuffers {
    SolverWorkspace solver;
    std::vector<double> outputs;
  };

  void initializeGridSearch(double targetDecoySizeRatio);
  void trainCpCnPair(CandidateCposCfrac& cpCnFold,
                     const options& pOptions,
                     const AlgIn& svmInput,
                     SvmTrainingBuffers& buffers);

  int mergeCpCnPairs(double selectionFdr,
                     const options& pOptions,
                     std::vector<std::vector<Scores> >& nestedTestScoresVec,
                     const vector<double>& cpos_vec,
                     const vector<double>& cfrac_vec);
  int doStep(const Normalizer* pNorm, double selectionFdr);
  // throws, or warns with --no-terminate, if a weight of a fold is not finite
  void checkWeightsFinite() const;
  // lists the input line and fold of every test PSM in testFolds_
  void recordTestFolds();

  void printSetWeights(ostream& weightStream, unsigned int set);
  void printRawSetWeights(ostream& weightStream,
                          unsigned int set,
                          const Normalizer* pNorm);

  void printAllWeightsColumns(ostream& weightStream);
  void printAllRawWeightsColumns(ostream& weightStream,
                                 const Normalizer* pNorm);
  void printAllWeightsColumns(
      std::vector<std::vector<double> > weightMatrix,
      ostream& weightStream);
};

#endif /*CROSSVALIDATION_H_*/
