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

#include "CrossValidation.h"
#include "PosteriorEstimator.h"
#include "Timer.h"

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
// checks cross validation convergence in case of quickValidation_
const double CrossValidation::requiredIncreaseOver2Iterations_ = 0.01;
// checks convergence in normal (not classic) training
const double CrossValidation::requiredIncreaseOver2IterationsDefault_ = 0.001;
// number of iterations in which Cpos and Cneg are selected; afterwards each
// fold keeps its pair (not in classic training)
const unsigned int CrossValidation::gridSearchIterations_ = 3u;

CrossValidation::CrossValidation(bool quickValidation,
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
                                 unsigned int numFolds)
    : quickValidation_(quickValidation),
      usePi0_(usePi0),
      reportPerformanceEachIteration_(reportPerformanceEachIteration),
      testFdr_(testFdr),
      selectionFdr_(selectionFdr),
      initialSelectionFdr_(initialSelectionFdr),
      selectedCpos_(selectedCpos),
      selectedCneg_(selectedCneg),
      niter_(niter),
      nestedXvalBins_(nestedXvalBins),
      trainBestPositive_(trainBestPositive),
      numThreads_(numThreads),
      skipNormalizeScores_(skipNormalizeScores),
      decoyFractionTraining_(decoyFractionTraining),
      numFolds_(numFolds) {
  // initialize selectionFdr_ and initialSelectionFdr_ if they have not been set
  // explicitly by the user.
  if (selectionFdr_ <= 0.0) {
    selectionFdr_ = testFdr_;
    if (initialSelectionFdr_ <= 0.0) {
      initialSelectionFdr_ = testFdr_;
    }
  }
}

/**
 * Sets up the SVM classifier:
 * - divide dataset into training and test sets for each fold
 * - set parameters (fdr, soft margin)
 *
 * Side effects: updates the following member variables:
 * - `svmInputs_`
 * - `trainScores_`, `testScores_` (in createXvalSetsBySpectrum)
 * - `weights_`
 * - `candidatesCpos_`, `candidatesCfrac_` and `classWeightsPerFold_` (in
 * `initializeGridSearch`)
 *
 * @param w_ vector of SVM weights
 * @return number of positives for initial setup
 */
int CrossValidation::preIterationSetup(Scores& fullset,
                                       SanityCheck& check,
                                       const Normalizer* pNorm) {
  assert(nestedXvalBins_ >= 1u);
  // One input set per (nested) cross validation bin, to be reused multiple
  // times
  const int numFeatures = static_cast<int>(FeatureNames::getNumFeatures());
  svmInputs_.assign(numFolds_ * nestedXvalBins_, AlgIn(numFeatures + 1));

  unsigned int decoysPerTarget = 1u;  // TODO: make this a command line argument
  {
    StageTimer stageTimer("cross-validation split");
    fullset.createXvalSetsBySpectrum(trainScores_, testScores_, numFolds_,
                                     decoyFractionTraining_, decoysPerTarget);
  }

  // initialize weights vector for all folds
  weights_ = vector<vector<double> >(
      numFolds_, vector<double>(FeatureNames::getNumFeatures() + 1));
  StageTimer stageTimer("initial direction");
  int numPositive =
      check.getInitDirection(testScores_, trainScores_, pNorm, weights_,
                               testFdr_, initialSelectionFdr_);

  initializeGridSearch(fullset.getTargetDecoySizeRatio());

  // the folds hold copies of all PSMs; postIterationProcessing merges them
  // back into fullset
  fullset.reset();
  return numPositive;
}

/**
 * Fill member variables with grid search parameters.
 *
 * Side effects: updates the candidatesCpos_, candidatesCfrac_ and
 * classWeightsPerFold_ member variables.
 *
 * @param targetDecoySizeRatio ratio of target to decoy PSMs
 */
void CrossValidation::initializeGridSearch(double targetDecoySizeRatio) {
  if (selectedCpos_ > 0) {
    candidatesCpos_.push_back(selectedCpos_);
  } else {
    candidatesCpos_.push_back(10);
    candidatesCpos_.push_back(1);
    candidatesCpos_.push_back(0.1);
    if (VERB > 0) {
      cerr << "Selecting Cpos by cross-validation." << endl;
    }
  }
  if (selectedCpos_ > 0 && selectedCneg_ > 0) {
    candidatesCfrac_.push_back(selectedCneg_ / selectedCpos_);
  } else {
    candidatesCfrac_.push_back(1.0 * targetDecoySizeRatio);
    candidatesCfrac_.push_back(3.0 * targetDecoySizeRatio);
    candidatesCfrac_.push_back(10.0 * targetDecoySizeRatio);
    if (VERB > 0) {
      cerr << "Selecting Cneg by cross-validation." << endl;
    }
  }

  // Quick validation trains only Cpos = Cneg = 1. The candidate lists are
  // replaced as well, since the selection after nested cross validation
  // loops over them and would otherwise pick a pair that was never trained.
  if (quickValidation_) {
    candidatesCpos_.assign(1, 1.0);
    candidatesCfrac_.assign(1, 1.0);
  }

  // Form cpos, cneg pairs per nested CV fold
  CandidateCposCfrac cpCnFold;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    for (int nestedSet = 0; nestedSet < nestedXvalBins_; nestedSet++) {
      for (double cpos : candidatesCpos_) {
        for (double cfrac : candidatesCfrac_) {
          cpCnFold.cpos = cpos;
          cpCnFold.cfrac = cfrac;
          cpCnFold.set = set;
          cpCnFold.nestedSet = nestedSet;
          cpCnFold.tp = 0;
          cpCnFold.ww.clear();
          for (int i = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
               i--;) {
            cpCnFold.ww.push_back(0);
          }
          classWeightsPerFold_.push_back(cpCnFold);
        }
      }
    }
  }
}

/**
 * Train the SVM using several cross validation iterations
 * @param pNorm Normalization object
 */
void CrossValidation::train(const Normalizer* pNorm) {
  if (VERB > 0) {
    cerr << "---Training with Cpos";
    if (selectedCpos_ > 0) {
      cerr << "=" << selectedCpos_;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", Cneg";
    if (selectedCneg_ > 0) {
      cerr << "=" << selectedCneg_;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", initial_fdr=" << initialSelectionFdr_;
    cerr << ", fdr=" << selectionFdr_ << endl;
  }

  // iterate
  int foundPositivesOldOld = 0, foundPositivesOld = 0, foundPositives = 0;
  for (unsigned int i = 0; i < niter_; i++) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << ":\t";
    }

    double selectionFdr = selectionFdr_;
    if (i == 0u) {
      selectionFdr = initialSelectionFdr_;
    }
    foundPositives = doStep(pNorm, selectionFdr);

    if (reportPerformanceEachIteration_) {
      int foundTestPositives = 0;
      for (size_t set = 0; set < numFolds_; ++set) {
        foundTestPositives +=
            testScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
      }
      if (VERB > 1) {
        std::cerr << "Found " << foundTestPositives << " test set PSMs with q<"
                  << testFdr_ << endl;
      }
    } else if (VERB > 1) {
      cerr << "Estimated " << foundPositives << " PSMs with q<" << testFdr_
           << endl;
    }

    if (VERB > 2) {
      printAllWeightsColumns(cerr);
    }
    if (VERB > 3) {
      printAllRawWeightsColumns(cerr, pNorm);
    }
    // stop when the estimated positives grew by less than a fraction over the
    // last two iterations; classic training runs all iterations unless
    // quick validation is used
    const double requiredIncrease = quickValidation_
                                        ? requiredIncreaseOver2Iterations_
                                        : requiredIncreaseOver2IterationsDefault_;
    if (foundPositives > 0 && foundPositivesOldOld > 0 &&
        (quickValidation_ || !classicTraining_) &&
        (static_cast<double>(foundPositives - foundPositivesOldOld) <=
         foundPositivesOldOld * requiredIncrease)) {
      if (VERB > 0) {
        std::cerr << "Performance increase over the last two iterations "
                  << "indicate that the algorithm has converged\n"
                  << "(" << foundPositives << " vs " << foundPositivesOldOld
                  << ")" << std::endl;
      }
      break;
    }
    foundPositivesOldOld = foundPositivesOld;
    foundPositivesOld = foundPositives;
  }
  if (VERB == 2) {
    printAllWeightsColumns(cerr);
  }
  if (VERB == 3) {
    printAllRawWeightsColumns(cerr, pNorm);
  }
  foundPositives = 0;
  for (size_t set = 0; set < numFolds_; ++set) {
    foundPositives +=
        testScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
  }
  if (VERB > 0) {
    std::cerr << "Found " << foundPositives << " test set PSMs with q<"
              << testFdr_ << "." << std::endl;
  }
}

/**
 * Executes a cross validation step
 * @param w_ list of the bins' normal vectors (in linear algebra sense) of the
 *        hyperplane from SVM
 * @return Estimation of number of true positives
 */
int CrossValidation::doStep(const Normalizer* pNorm, double selectionFdr) {
  // Setup
  options pOptions;
  pOptions.lambda = 1.0;
  pOptions.lambda_u = 1.0;
  pOptions.epsilon = EPSILON;
  pOptions.cgitermax = CGITERMAX;
  pOptions.mfnitermax = MFNITERMAX;
  pOptions.directSolve = !classicTraining_;
  int estTruePos = 0;

  // for determining an appropriate positive training set, the decoys+1 in the
  // FDR estimates is too restrictive for small datasets
  bool skipDecoysPlusOne = true;
  {
    StageTimer stageTimer("rescore training sets");
    for (std::size_t set = 0; set < numFolds_; ++set) {
      if (trainScoresCurrent_) {
        // already scored and sorted with weights_ at the end of the last step
        trainScores_[set].calcQvals(selectionFdr, skipDecoysPlusOne);
      } else {
        trainScores_[set].calcScoresAndQvals(weights_[set], selectionFdr,
                                             skipDecoysPlusOne);
      }
    }
  }

  // Below implements the series of speedups detailed in the following:
  // ////////////////////////////////
  // Speeding Up Percolator
  // John T. Halloran, Hantian Zhang, Kaan Kara, Cédric Renggli, Matthew The, Ce
  // Zhang, David M. Rocke, Lukas Käll, and William Stafford Noble
  //   Journal of Proteome Research 2019 18 (9), 3353-3359
  // ////////////////////////////////
  // Note that the implementation further improves on the speedups in the paper
  // by:
  //   -implementing a single threadpool using OMP for SVM training per each
  //   cpos,cneg pair per nested CV fold -has a much smaller memory footprint by
  //   fixing memory leaks in L2_SVM_MFN and more efficient validation of the
  //   learned SVM parameters
  // ////

  // Create SVM input data for parallelization
  std::vector<std::vector<Scores> > nestedTestScoresVec;
  {
    StageTimer stageTimer("build SVM input");
    for (std::size_t set = 0; set < numFolds_; ++set) {
      // Without nested cross validation, the training set is also used for
      // validation (sub-optimal cross validation).
      std::vector<Scores> nestedTrainScores, nestedTestScores;
      if (nestedXvalBins_ > 1) {
        trainScores_[set].createXvalSetsBySpectrum(
            nestedTrainScores, nestedTestScores, nestedXvalBins_);
      }
      nestedTestScoresVec.push_back(std::move(nestedTestScores));
      // Set SVM input data for L2-SVM-MFN
      for (std::size_t nestedFold = 0; nestedFold < nestedXvalBins_;
           ++nestedFold) {
        AlgIn& svmInput = svmInputs_[set * nestedXvalBins_ + nestedFold];
        Scores& nestedTrainSet = nestedTrainScores.empty()
                                     ? trainScores_[set]
                                     : nestedTrainScores[nestedFold];
        if ((VERB > 2) && (nestedFold == 0)) {
          cerr << "Split " << set + 1 << ": Training with "
               << svmInput.positives << " positives and "
               << svmInput.negatives << " negatives" << std::endl;
        }
        nestedTrainSet.generateNegativeTrainingSet(svmInput, 1.0);
        nestedTrainSet.generatePositiveTrainingSet(
            svmInput, selectionFdr, 1.0, trainBestPositive_);
        if ((VERB > 2) && (nestedXvalBins_ > 1u)) {
          cerr << "Split " << set + 1 << " nested fold " << nestedFold
               << ": Training with " << svmInput.positives << " positives and "
               << svmInput.negatives << " negatives" << std::endl;
        }
      }
    }
  }

  {
    StageTimer stageTimer("SVM training");
#pragma omp parallel for schedule(dynamic, 1) ordered
    for (int pairIdx = 0; pairIdx < classWeightsPerFold_.size(); pairIdx++) {
      CandidateCposCfrac* cpCnFold = &classWeightsPerFold_[pairIdx];
      if (!isSelectedPair(*cpCnFold)) {
        continue;
      }
      const AlgIn& svmInput =
          svmInputs_[cpCnFold->set * nestedXvalBins_ +
                     static_cast<unsigned int>(cpCnFold->nestedSet)];
      trainCpCnPair(*cpCnFold, pOptions, svmInput);
    }
  }

  StageTimer stageTimer("grid evaluation");
  estTruePos = mergeCpCnPairs(selectionFdr, pOptions, nestedTestScoresVec,
                              candidatesCpos_, candidatesCfrac_);
  return estTruePos;
}

/**
 * Train SVM over a single (cpos, cneg) pair
 * @param cpCnFold contains cpos, cneg pair and SVM learned weights
 * @param pOptions options for the SVM algorithm
 * @param svmInput training data for this particular nested CV fold
 */
void CrossValidation::trainCpCnPair(CandidateCposCfrac& cpCnFold,
                                    const options& pOptions,
                                    const AlgIn& svmInput) {
  double cpos = cpCnFold.cpos;
  double cfrac = cpCnFold.cfrac;

  if (VERB > 3)
    cerr << "- cross-validation with Cpos=" << cpos << ", Cneg=" << cfrac * cpos
         << endl;
  // The SVM algorithm (see ssl.cpp) starts from the weights this pair had
  // after the previous iteration (zero in the first one), which are close
  // to the new solution, or from zero with classic training.
  std::vector<double>& weights = cpCnFold.ww;
  std::vector<double> outputs(
      static_cast<std::size_t>(svmInput.positives + svmInput.negatives), 0.0);
  if (classicTraining_) {
    std::fill(weights.begin(), weights.end(), 0.0);
  } else {
    for (int i = 0; i < svmInput.m; ++i) {
      const double* x = svmInput.row(i);
      double o = 0.0;
      for (int k = 0; k < svmInput.n; ++k) o += x[k] * weights[k];
      outputs[i] = o;
    }
  }
  L2_SVM_MFN(svmInput, pOptions, weights, outputs, cpos, cfrac * cpos);
}

/**
 * Validate and merge weights learned per cpos,cneg pairs per nested CV fold per
 * CV fold
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
 * @param nestedTestScoresVec 2D vector, test sets per nested CV fold per CV
 * fold
 */
int CrossValidation::mergeCpCnPairs(
    double selectionFdr,
    const options& pOptions,
    vector<vector<Scores> >& nestedTestScoresVec,
    const vector<double>& cposCandidates,
    const vector<double>& cfracCandidates) {
  // for determining the number of positives, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true;

  vector<int> bestTruePoses(numFolds_, 0);
  vector<double> bestCposes(numFolds_, 1);
  vector<double> bestCfracs(numFolds_, 1);
  
  // Validate learned parameters per (cpos,cneg) pair per nested CV fold
  // Note: this cannot be done in trainCpCnPair without setting a critical
  // pragma, due to the
  //       scoring calculation in calcScores.
  unsigned int numCpCnPairsPerSet =
      static_cast<unsigned int>(classWeightsPerFold_.size() / numFolds_);

  // The validation counts of all pairs are independent, so they are computed
  // in parallel. Without nested cross validation the training sets are the
  // validation sets.
  auto testSet = [&](unsigned int set,
                     unsigned int nestedSet) -> const Scores& {
    const std::vector<Scores>& nestedTestScores = nestedTestScoresVec[set];
    return nestedTestScores.empty() ? trainScores_[set]
                                    : nestedTestScores[nestedSet];
  };
  const int numTestSets = static_cast<int>(numFolds_ * nestedXvalBins_);
  std::vector<Scores::FeatureRows> testRows(numTestSets);
#pragma omp parallel for schedule(dynamic, 1)
  for (int idx = 0; idx < numTestSets; ++idx) {
    testSet(idx / nestedXvalBins_, idx % nestedXvalBins_)
        .getFeatureRows(testRows[idx]);
  }

  std::vector<int> truePositives(classWeightsPerFold_.size());
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
#pragma omp parallel for schedule(dynamic, 1)
  for (int pairIdx = 0; pairIdx < static_cast<int>(classWeightsPerFold_.size());
       ++pairIdx) {
    const CandidateCposCfrac& cpCnPair = classWeightsPerFold_[pairIdx];
    const unsigned int nestedSet =
        static_cast<unsigned int>(cpCnPair.nestedSet);
    std::vector<std::pair<double, bool> > combined;
    if (!isSelectedPair(cpCnPair)) {
      truePositives[pairIdx] = -1;  // not trained in this iteration
      continue;
    }
    truePositives[pairIdx] =
        testSet(cpCnPair.set, nestedSet)
            .countTargetsBelowFdr(
                testRows[cpCnPair.set * nestedXvalBins_ + nestedSet],
                cpCnPair.ww, testFdr_, skipDecoysPlusOne, combined);
  }

  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    unsigned int a = set * numCpCnPairsPerSet;
    unsigned int b = (set + 1) * numCpCnPairsPerSet;
    int tp = 0;
    std::vector<CandidateCposCfrac>::iterator itCpCnPair;
    std::map<std::pair<double, double>, int> intermediateResults;
    for (itCpCnPair = classWeightsPerFold_.begin() + a;
         itCpCnPair < classWeightsPerFold_.begin() + b; itCpCnPair++) {
      tp = truePositives[itCpCnPair - classWeightsPerFold_.begin()];
      intermediateResults[std::pair<double,double>{itCpCnPair->cpos,
                                                   itCpCnPair->cfrac}] += tp;
      itCpCnPair->tp = tp;
      if (nestedXvalBins_ <= 1) {
        if (tp >= bestTruePoses[set]) {
          bestTruePoses[set] = tp;
          weights_[set] = itCpCnPair->ww;
          bestCposes[set] = itCpCnPair->cpos;
          bestCfracs[set] = itCpCnPair->cfrac;
        }
      }
    }
    if (nestedXvalBins_ > 1) {  // Check nestedXvalBins, which collapse
                                // (accumulate) tp estimated for each CV bin
      // Now check which achieved best performance among cpos, cneg pairs
      std::vector<double>::const_iterator itCpos = cposCandidates.begin();
      for (; itCpos != cposCandidates.end(); ++itCpos) {
        double cpos = *itCpos;
        std::vector<double>::const_iterator itCfrac = cfracCandidates.begin();
        for (; itCfrac != cfracCandidates.end(); ++itCfrac) {
          double cfrac = *itCfrac;
          tp = intermediateResults[std::pair<double,double>{cpos, cfrac}];
          if (tp >= bestTruePoses[set]) {
            bestTruePoses[set] = tp;
            bestCposes[set] = cpos;
            bestCfracs[set] = cfrac;
          }
        }
      }
    }
  }

  if (nestedXvalBins_ > 1) {
#pragma omp parallel for schedule(dynamic, 1) ordered
    for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
      AlgIn& svmInput = svmInputs_[set * nestedXvalBins_];
      trainScores_[set].generateNegativeTrainingSet(svmInput, 1.0);
      trainScores_[set].generatePositiveTrainingSet(svmInput, selectionFdr,
                                                    1.0, trainBestPositive_);

      // the SVM algorithm (see ssl.cpp) starts from zero weights and outputs
      std::fill(weights_[set].begin(), weights_[set].end(), 0.0);
      std::vector<double> outputs(
          static_cast<std::size_t>(svmInput.positives + svmInput.negatives),
          0.0);
      L2_SVM_MFN(svmInput, pOptions, weights_[set], outputs, bestCposes[set],
                 bestCposes[set] * bestCfracs[set]);
    }
  }

  // After the first iterations the selected pairs hardly change, and their
  // number of positives differs from the other pairs' only by noise. From
  // then on only the selected pair of each fold is trained.
  ++iteration_;
  if (!classicTraining_ && !gridFixed_ && iteration_ >= gridSearchIterations_ &&
      classWeightsPerFold_.size() > numFolds_ * nestedXvalBins_) {
    gridFixed_ = true;
    selectedPairs_.resize(numFolds_);
    for (unsigned int set = 0; set < numFolds_; ++set) {
      selectedPairs_[set] = std::make_pair(bestCposes[set], bestCfracs[set]);
      if (VERB > 1) {
        std::cerr << "Split " << set + 1 << ": using Cpos=" << bestCposes[set]
                  << ", Cneg=" << bestCposes[set] * bestCfracs[set]
                  << " in the remaining iterations" << std::endl;
      }
    }
  }

  double bestTruePos = 0;
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    bestTruePos +=
        trainScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
  }
  trainScoresCurrent_ = true;
  // each PSM is in numFolds_ - 1 training sets, or in the only one
  return static_cast<int>(bestTruePos / std::max(1u, numFolds_ - 1));
}

void CrossValidation::postIterationProcessing(Scores& fullset,
                                              SanityCheck& check) {
  if (decoyFractionTraining_ < 1.0) {
    // with RESET 1 decoy: decoysPerTarget = 1 and decoyFractionTraining = 0.5,
    // then factor = 1.5 and testNullTargetWinProb = 2/3 (and decoyFactor = 2)
    //
    // with RESET 2 decoys: decoysPerTarget = 2 and decoyFractionTraining = 0.5,
    // then factor = 2 and testNullTargetWinProb = 1/2 (and decoyFactor = 1)
    unsigned int decoysPerTarget =
        1u;  // TODO: make this a command line argument
    double factor =
        (decoysPerTarget + 1) - decoyFractionTraining_ * decoysPerTarget;
    double testNullTargetWinProb = 1.0 / factor;

    for (std::size_t set = 0; set < numFolds_; ++set) {
      Scores newTestScores(false);
      for (ScoreHolder& sh : testScores_[set]) {
        if (sh.label == LabelType::DECOY)
          continue;

        if (sh.label == LabelType::PSEUDO_TARGET) {
          sh.label = LabelType::DECOY;
        }
        newTestScores.addScoreHolder(sh);
      }
      newTestScores.setNullTargetWinProb(testNullTargetWinProb);
      newTestScores.recalculateSizes();
      testScores_[set] = std::move(newTestScores);
    }
    fullset.setNullTargetWinProb(testNullTargetWinProb);
  }

  if (!check.validateDirection(testScores_, trainScores_, weights_)) {
    for (std::size_t set = 0; set < numFolds_; ++set) {
      testScores_[set].calcScoresAndQvals(weights_[0], selectionFdr_);
    }
  }
  fullset.merge(testScores_, selectionFdr_, skipNormalizeScores_, weights_);
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printSetWeights(ostream& weightStream, unsigned int set) {
  weightStream << DataSet::getFeatureNames().getFeatureNames() << "\tm0"
               << std::endl;
  weightStream.precision(3);
  weightStream << weights_[set][0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << weights_[set][ix];
  }
  weightStream << endl;
}

/**
 * Prints the weights of the raw vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printRawSetWeights(ostream& weightStream,
                                         unsigned int set,
                                         const Normalizer* pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(weights_[set], ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

// used to save weights for reuse as weight input file
void CrossValidation::printAllWeights(ostream& weightStream,
                                      const Normalizer* pNorm) {
  weightStream << "# This file contains the weights from each cross validation "
                  "bin from percolator training"
               << endl;
  weightStream << "# First line is the feature names, followed by normalized "
                  "weights, and the raw weights of bin 1"
               << endl;
  weightStream << "# This is repeated for the other bins" << endl;
  for (unsigned int ix = 0; ix < numFolds_; ++ix) {
    printSetWeights(weightStream, ix);
    printRawSetWeights(weightStream, ix, pNorm);
  }
}

void CrossValidation::getAvgWeights(std::vector<double>& weights,
                                    const Normalizer* pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);

  weights.resize(FeatureNames::getNumFeatures() + 1);
  for (size_t set = 0; set < weights_.size(); ++set) {
    pNorm->unnormalizeweight(weights_[set], ww);
    for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
      weights[ix] += ww[ix] / static_cast<double>(weights_.size());
    }
  }
}

void CrossValidation::printAllWeightsColumns(
    std::vector<std::vector<double> > weightMatrix,
    ostream& outputStream) {
  // write to intermediate stream to prevent the fixed precision from sticking
  ostringstream weightStream;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    weightStream << " Split" << set + 1 << '\t';  // right-align with weights
  }
  weightStream << "FeatureName" << std::endl;
  size_t numRows = FeatureNames::getNumFeatures() + 1;
  for (unsigned int ix = 0; ix < numRows; ix++) {
    for (unsigned int set = 0; set < numFolds_; ++set) {
      // align positive and negative weights
      if (weightMatrix[set][ix] >= 0)
        weightStream << " ";
      weightStream << fixed << setprecision(4) << weightMatrix[set][ix];
      weightStream << "\t";
    }
    if (ix == numRows - 1) {
      weightStream << "m0";
    } else {
      weightStream << DataSet::getFeatureNames().getFeatureName(ix);
    }
    weightStream << std::endl;
  }
  outputStream << weightStream.str();
}

void CrossValidation::printAllWeightsColumns(ostream& weightStream) {
  std::cerr << "Learned normalized SVM weights for the " << numFolds_
            << " cross-validation splits:" << std::endl;
  printAllWeightsColumns(weights_, weightStream);
}

void CrossValidation::printAllRawWeightsColumns(ostream& weightStream,
                                                const Normalizer* pNorm) {
  std::cerr << "Learned raw SVM weights for the " << numFolds_
            << " cross-validation splits:" << std::endl;
  std::vector<std::vector<double> > ww = weights_;
  for (size_t set = 0; set < weights_.size(); ++set) {
    pNorm->unnormalizeweight(weights_[set], ww[set]);
  }
  printAllWeightsColumns(ww, weightStream);
}
