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

// number of folds for cross validation
const unsigned int CrossValidation::numFolds_ = 3u;
#ifdef _OPENMP
const unsigned int CrossValidation::numAlgInObjects_ = CrossValidation::numFolds_;
#else
const unsigned int CrossValidation::numAlgInObjects_ = 1u;
#endif
// checks cross validation convergence in case of quickValidation_
const double CrossValidation::requiredIncreaseOver2Iterations_ = 0.01; 

CrossValidation::CrossValidation(bool quickValidation, 
  bool reportPerformanceEachIteration, double testFdr, double selectionFdr, 
  double selectedCpos, double selectedCneg, int niter, bool usePi0) :
    quickValidation_(quickValidation), usePi0_(usePi0),
    reportPerformanceEachIteration_(reportPerformanceEachIteration), 
    testFdr_(testFdr), selectionFdr_(selectionFdr), 
    selectedCpos_(selectedCpos), selectedCneg_(selectedCneg), niter_(niter) {}


CrossValidation::~CrossValidation() { 
  for (unsigned int set = 0; set < numAlgInObjects_; ++set) {
    if (svmInputs_[set]) {
      delete svmInputs_[set];
    }
    svmInputs_[set] = NULL;
  }
}

/** 
 * Sets up the SVM classifier: 
 * - divide dataset into training and test sets for each fold
 * - set parameters (fdr, soft margin)
 * @param w_ vector of SVM weights
 * @return number of positives for initial setup
 */
int CrossValidation::preIterationSetup(Scores& fullset, SanityCheck* pCheck, 
                                       Normalizer* pNorm) {
  // initialize weights vector for all folds
  w_ = vector<vector<double> >(numFolds_, 
           vector<double> (FeatureNames::getNumFeatures() + 1));
  
  // One input set, to be reused multiple times
  for (unsigned int set = 0; set < numAlgInObjects_; ++set) {
    svmInputs_.push_back(new AlgIn(fullset.size(), FeatureNames::getNumFeatures() + 1));
    assert( svmInputs_.back() );
  }
  
  if (selectedCpos_ >= 0 && selectedCneg_ >= 0) {
    trainScores_.resize(numFolds_, Scores(usePi0_));
    testScores_.resize(numFolds_, Scores(usePi0_));
    
    fullset.createXvalSetsBySpectrum(trainScores_, testScores_, numFolds_);
    
    if (selectionFdr_ <= 0.0) {
      selectionFdr_ = testFdr_;
    }
    if (selectedCpos_ > 0) {
      candidatesCpos_.push_back(selectedCpos_);
    } else {
      candidatesCpos_.push_back(10);
      candidatesCpos_.push_back(1);
      candidatesCpos_.push_back(0.1);
      if (VERB > 0) {
        cerr << "selecting cpos by cross validation" << endl;
      }
    }
    if (selectedCpos_ > 0 && selectedCneg_ > 0) {
      candidatesCfrac_.push_back(selectedCneg_ / selectedCpos_);
    } else {
      candidatesCfrac_.push_back(1.0 * fullset.getTargetDecoySizeRatio());
      candidatesCfrac_.push_back(3.0 * fullset.getTargetDecoySizeRatio());
      candidatesCfrac_.push_back(10.0 * fullset.getTargetDecoySizeRatio());
      if (VERB > 0) {
        cerr << "selecting cneg by cross validation" << endl;
      }
    }
    return pCheck->getInitDirection(testScores_, trainScores_, pNorm, w_, 
                                    testFdr_);
  } else {
    vector<Scores> myset(1, fullset);
    cerr << "B" << endl;
    return pCheck->getInitDirection(myset, myset, pNorm, w_, testFdr_);
  }
}

/** 
 * Train the SVM using several cross validation iterations
 * @param pNorm Normalization object
 */
void CrossValidation::train(Normalizer * pNorm) {

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
    cerr << ", fdr=" << selectionFdr_ << endl;
  }
  
  // iterate
  int foundPositivesOldOld = 0, foundPositivesOld = 0, foundPositives = 0; 
  for (unsigned int i = 0; i < niter_; i++) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << " :\t";
    }
    
    bool updateDOC = true;
    foundPositives = doStep(updateDOC);
    
    if (VERB > 1) {
      cerr << "After the iteration step, " << foundPositives
          << " target PSMs with q<" << testFdr_
          << " were estimated by cross validation" << endl;
    }
    if (VERB > 2) {
      cerr << "Obtained weights" << endl;
      printAllWeights(cerr, pNorm);
    }
    if (foundPositives > 0 && foundPositivesOldOld > 0 && quickValidation_ &&
           (static_cast<double>(foundPositives - foundPositivesOldOld) <= 
           foundPositivesOldOld * requiredIncreaseOver2Iterations_)) {
      if (VERB > 1) {
        std::cerr << "Performance increase over the last two iterations " <<
            "indicate that the algorithm has converged\n" <<
            "(" << foundPositives << " vs " << foundPositivesOldOld << ")" << 
            std::endl;
      }
      break;
    }
    foundPositivesOldOld = foundPositivesOld;    
    foundPositivesOld = foundPositives;
  }
  if (VERB == 2) {
    cerr
    << "Obtained weights (only showing weights of first cross validation set)"
    << endl;
    printSetWeights(cerr, 0, pNorm);
  }
  foundPositives = 0;
  for (size_t set = 0; set < numFolds_; ++set) {
    if (DataSet::getCalcDoc()) {
      testScores_[set].getDOC().copyDOCparameters(trainScores_[set].getDOC());
      testScores_[set].setDOCFeatures();
    }
    foundPositives += testScores_[set].calcScores(w_[set], testFdr_);
  }
  if (VERB > 0) {
    std::cerr << "After all training done, " << foundPositives << 
                 " target PSMs with q<" << testFdr_ << 
                 " were found when measuring on the test set" << std::endl;
  }  
}


/** 
 * Executes a cross validation step
 * @param w_ list of the bins' normal vectors (in linear algebra sense) of the 
 *        hyperplane from SVM
 * @param updateDOC boolean deciding to recalculate retention features 
 *        @see DescriptionOfCorrect
 * @return Estimation of number of true positives
 */
int CrossValidation::doStep(bool updateDOC) {
  // Setup
  struct options* pOptions = new options;
  pOptions->lambda = 1.0;
  pOptions->lambda_u = 1.0;
  pOptions->epsilon = EPSILON;
  pOptions->cgitermax = CGITERMAX;
  pOptions->mfnitermax = MFNITERMAX;
  int estTruePos = 0;
  if (!quickValidation_) {
  #pragma omp parallel for schedule(dynamic, 1)
    for (int set = 0; set < numFolds_; ++set) {
      struct vector_double* pWeights = new vector_double;
      pWeights->d = FeatureNames::getNumFeatures() + 1;
      pWeights->vec = new double[pWeights->d];
      
      double bestCpos = 1, bestCfrac = 1;
      
      int estTruePosFold = processSingleFold(set, updateDOC, candidatesCpos_, 
                                      candidatesCfrac_, bestCpos, bestCfrac, 
                                      pWeights, pOptions);
      #pragma omp critical (add_tps)
      {
        estTruePos += estTruePosFold;
      }
      
      delete[] pWeights->vec;
      delete pWeights;
    }
  } else {
    struct vector_double* pWeights = new vector_double;
    pWeights->d = FeatureNames::getNumFeatures() + 1;
    pWeights->vec = new double[pWeights->d];
    
    double bestCpos = 1, bestCfrac = 1;
    
    // Use limited internal cross validation, i.e take the cpos and cfrac 
    // values of the first bin and use it for the subsequent bins 
    estTruePos += processSingleFold(0, updateDOC, candidatesCpos_, 
                                    candidatesCfrac_, bestCpos, bestCfrac, 
                                    pWeights, pOptions);
    vector<double> cp(1, bestCpos), cf(1, bestCfrac);
  #pragma omp parallel for schedule(dynamic, 1)
    for (int set = 1; set < numFolds_; ++set) {
      int estTruePosFold = processSingleFold(set, updateDOC, cp, cf, bestCpos, 
                                      bestCfrac, pWeights, pOptions);
      #pragma omp critical (add_tps)
      {
        estTruePos += estTruePosFold;
      }  
    }
    delete[] pWeights->vec;
    delete pWeights;
  }
  delete pOptions;
  return estTruePos / (numFolds_ - 1);
}

/** 
 * Train one of the crossvalidation bins 
 * @param set identification number of the bin that is processed
 * @param updateDOC boolean deciding to calculate retention features 
 *        @see DescriptionOfCorrect
 * @param cposCandidates candidate soft margin parameters for positives
 * @param cfracCandidates candidate soft margin parameters for fraction neg/pos
 * @param bestCpos best soft margin parameter for positives
 * @param bestCfrac best soft margin parameter for fraction neg/pos
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
*/
int CrossValidation::processSingleFold(unsigned int set, bool updateDOC, 
    const vector<double>& cposCandidates, const vector<double>& cfracCandidates, 
    double &bestCpos, double &bestCfrac, vector_double* pWeights, 
    options * pOptions) {
  int bestTruePos = 0;
  if (VERB > 2) {
    cerr << "cross validation - fold " << set + 1 << " out of "
         << numFolds_ << endl;
  }
  
  vector<double> ww = w_[set]; // SVM weights initial guess and result holder
  vector<double> bestW = w_[set]; // SVM weights with highest true pos estimate
  trainScores_[set].calcScores(ww, selectionFdr_);
  if (DataSet::getCalcDoc() && updateDOC) {
    trainScores_[set].recalculateDescriptionOfCorrect(selectionFdr_);
  }
  
  AlgIn* svmInput = svmInputs_[set % numAlgInObjects_];
  
  trainScores_[set].generateNegativeTrainingSet(*svmInput, 1.0);
  trainScores_[set].generatePositiveTrainingSet(*svmInput, selectionFdr_, 1.0);
  if (VERB > 2) {
    cerr << "Calling with " << svmInput->positives << " positives and "
         << svmInput->negatives << " negatives\n";
  }
  
  // Create storage vector for SVM algorithm
  struct vector_double* Outputs = new vector_double;
  size_t numInputs = svmInput->positives + svmInput->negatives;
  Outputs->vec = new double[numInputs];
  Outputs->d = numInputs;
  
  // Find soft margin parameters with highest estimate of true positives
  std::vector<double>::const_iterator itCpos = cposCandidates.begin();
  for ( ; itCpos != cposCandidates.end(); ++itCpos) {
    double cpos = *itCpos;  
    std::vector<double>::const_iterator itCfrac = cfracCandidates.begin();
    for ( ; itCfrac != cfracCandidates.end(); ++itCfrac) {
      double cfrac = *itCfrac;
      if (VERB > 2) cerr << "-cross validation with cpos=" << cpos
          << ", cfrac=" << cfrac << endl;
      int tp = 0;
      for (int ix = 0; ix < pWeights->d; ix++) {
        pWeights->vec[ix] = 0;
      }
      for (int ix = 0; ix < Outputs->d; ix++) {
        Outputs->vec[ix] = 0;
      }
      svmInput->setCost(cpos, cpos * cfrac);
      
      // Call SVM algorithm (see ssl.cpp)
      L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs);
      
      for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
        ww[i] = pWeights->vec[i];
      }
      // sub-optimal cross validation (better would be to measure
      // performance on a set disjoint of the training set)
      tp = trainScores_[set].calcScores(ww, testFdr_);
      if (VERB > 2) {
        cerr << "- cross validation estimates " << tp
             << " target PSMs over " << testFdr_ * 100 << "% FDR level"
             << endl;
      }
      if (tp >= bestTruePos) {
        if (VERB > 2) {
          cerr << "Better than previous result, store this" << endl;
        }
        bestTruePos = tp;
        bestW = ww;
        bestCpos = cpos;
        bestCfrac = cfrac;
      }
    }
    if (VERB > 2) {
      std::cerr << "cross validation estimates " << 
          bestTruePos / (numFolds_-1) << " target PSMs with q<" << testFdr_ <<
          " for hyperparameters Cpos=" << bestCpos << 
          ", Cneg=" << bestCfrac * bestCpos << std::endl;
    }
  }
  w_[set] = bestW;
  delete[] Outputs->vec;
  delete Outputs;
  return bestTruePos;
}

void CrossValidation::postIterationProcessing(Scores& fullset,
                                              SanityCheck* pCheck) {
  if (!pCheck->validateDirection(w_)) {
    fullset.calcScores(w_[0], selectionFdr_);
  }
  if (VERB > 0) {
    std::cerr << "Merging results from " << testScores_.size() << 
                 " datasets" << std::endl;
  }
  fullset.merge(testScores_, selectionFdr_);
  for (unsigned int i = 0; i < numFolds_; ++i) {
    testScores_[i].deleteContiguousMemoryBlock();
  }
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printSetWeights(ostream & weightStream, unsigned int set, 
                                      Normalizer * pNorm) {
  weightStream << "# first line contains normalized weights, " <<
                  "second line the raw weights" << std::endl;
  weightStream << DataSet::getFeatureNames().getFeatureNames() << 
                  "\tm0" << std::endl;
  weightStream.precision(3);
  weightStream << w_[set][0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << w_[set][ix];
  }
  weightStream << endl;
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(w_[set], ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

void CrossValidation::getAvgWeights(std::vector<double>& weights, 
                                    Normalizer * pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  
  weights.resize(FeatureNames::getNumFeatures() + 1);
  for (size_t set = 0; set < w_.size(); ++set) {
    pNorm->unnormalizeweight(w_[set], ww);
    for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
      weights[ix] += ww[ix] / w_.size();
    }
  }
}

void CrossValidation::printAllWeights(ostream & weightStream, 
                                      Normalizer * pNorm) {
  for (unsigned int ix = 0; ix < numFolds_; ++ix) {
    printSetWeights(weightStream, ix, pNorm);
  }
}

void CrossValidation::printDOC() {
  cerr << "For the cross validation sets the average deltaMass are ";
  for (size_t ix = 0; ix < testScores_.size(); ix++) {
    cerr << testScores_[ix].getDOC().getAvgDeltaMass() << " ";
  }
  cerr << "and average pI are ";
  for (size_t ix = 0; ix < testScores_.size(); ix++) {
    cerr << testScores_[ix].getDOC().getAvgPI() << " ";
  }
  cerr << endl;
}
