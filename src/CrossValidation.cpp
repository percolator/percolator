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
  double initialSelectionFdr, double selectedCpos, double selectedCneg, int niter, bool usePi0,
  int nestedXvalBins, bool trainBestPositive) :
    quickValidation_(quickValidation), usePi0_(usePi0),
    reportPerformanceEachIteration_(reportPerformanceEachIteration), 
    testFdr_(testFdr), selectionFdr_(selectionFdr), initialSelectionFdr_(initialSelectionFdr),
    selectedCpos_(selectedCpos), selectedCneg_(selectedCneg), niter_(niter),
    nestedXvalBins_(nestedXvalBins), trainBestPositive_(trainBestPositive) {}


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
    Normalizer* pNorm, FeatureMemoryPool& featurePool) {
  assert(nestedXvalBins_ >= 1u);
  
  // initialize weights vector for all folds
  w_ = vector<vector<double> >(numFolds_, 
           vector<double> (FeatureNames::getNumFeatures() + 1));
  
  // One input set, to be reused multiple times
  for (unsigned int set = 0; set < numAlgInObjects_; ++set) {
    svmInputs_.push_back(new AlgIn(fullset.size(), FeatureNames::getNumFeatures() + 1));
    assert( svmInputs_.back() );
  }
  
  int numPositive = 0;
  if (selectedCpos_ >= 0 && selectedCneg_ >= 0) {
    trainScores_.resize(numFolds_, Scores(usePi0_));
    testScores_.resize(numFolds_, Scores(usePi0_));
    
    fullset.createXvalSetsBySpectrum(trainScores_, testScores_, numFolds_, featurePool);
    
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
        cerr << "Selecting Cpos by cross-validation." << endl;
      }
    }
    if (selectedCpos_ > 0 && selectedCneg_ > 0) {
      candidatesCfrac_.push_back(selectedCneg_ / selectedCpos_);
    } else {
      candidatesCfrac_.push_back(1.0 * fullset.getTargetDecoySizeRatio());
      candidatesCfrac_.push_back(3.0 * fullset.getTargetDecoySizeRatio());
      candidatesCfrac_.push_back(10.0 * fullset.getTargetDecoySizeRatio());
      if (VERB > 0) {
        cerr << "Selecting Cneg by cross-validation." << endl;
      }
    }
    numPositive = pCheck->getInitDirection(testScores_, trainScores_, pNorm, w_, 
                                           testFdr_, initialSelectionFdr_);
  } else {
    vector<Scores> myset(1, fullset);
    numPositive = pCheck->getInitDirection(myset, myset, pNorm, w_, testFdr_, initialSelectionFdr_);
  }

  // Form cpos, cneg, set pairs
  cpCnTriple cpCnFold;
  for (int set = 0; set < numFolds_; ++set) {
    for (int nestedSet = 0; nestedSet < nestedXvalBins_; nestedSet++) {
      if(!quickValidation_) {
      std::vector<double>::const_iterator itCpos = candidatesCpos_.begin();
      for ( ; itCpos != candidatesCpos_.end(); ++itCpos) {
	double cpos = *itCpos;
	std::vector<double>::const_iterator itCfrac = candidatesCfrac_.begin();
	for ( ; itCfrac != candidatesCfrac_.end(); ++itCfrac) {
	  double cfrac = *itCfrac;
	  cpCnFold.cpos = cpos;
	  cpCnFold.cfrac = cfrac;
	  cpCnFold.set = set;
	  cpCnFold.nestedSet = nestedSet;
	  cpCnFold.tp = 0;
	  for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
	    cpCnFold.ww.push_back(0);
	  }	  
	  classWeightsPerFold_.push_back(cpCnFold);
	}
      }
      } else {
	  cpCnFold.cpos = 1;
	  cpCnFold.cfrac = 1;
	  cpCnFold.set = set;
	  cpCnFold.nestedSet = nestedSet;
	  cpCnFold.tp = 0;
	  for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
	    cpCnFold.ww.push_back(0);
	  }	  
	  classWeightsPerFold_.push_back(cpCnFold);
	}
    }
  }
  // cout << "Num (cp,cfrac,set,nestedFold) = " << classWeightsPerFold_.size() << "\n";
  // for (int trip = 0; trip < classWeightsPerFold_.size(); trip++){
  //   cout << classWeightsPerFold_[trip].cpos << ", " << classWeightsPerFold_[trip].cfrac  << ", " << classWeightsPerFold_[trip].set << ", " << classWeightsPerFold_[trip].nestedSet << "\n";
  // }
  
  if (DataSet::getCalcDoc()) {
    for (int set = 0; set < numFolds_; ++set) {
      trainScores_[set].calcScores(w_[set], selectionFdr_);
    }
  #pragma omp parallel for schedule(dynamic, 1)
    for (int set = 0; set < numFolds_; ++set) {
      trainScores_[set].recalculateDescriptionOfCorrect(selectionFdr_);
      testScores_[set].getDOC().copyDOCparameters(trainScores_[set].getDOC());
      testScores_[set].setDOCFeatures(pNorm);
    }
  }
  
  return numPositive;
}

/** 
 * Train the SVM using several cross validation iterations
 * @param pNorm Normalization object
 */
void CrossValidation::train(Normalizer* pNorm) {

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
    
    bool updateDOC = true;
    double selectionFdr = selectionFdr_;
    if (i == 0u) {
      selectionFdr = initialSelectionFdr_;
    }
    foundPositives = doStep(updateDOC, pNorm, selectionFdr);
    
    if (reportPerformanceEachIteration_) {
      int foundTestPositives = 0;
      for (size_t set = 0; set < numFolds_; ++set) {
        foundTestPositives += testScores_[set].calcScores(w_[set], testFdr_);
      }
      if (VERB > 1) {
        cerr << "Found " << foundTestPositives << " test set PSMs with q<" 
             << testFdr_ << endl;
      }
    } else if (VERB > 1) {
      cerr << "Estimated " << foundPositives << " PSMs with q<" << testFdr_ << endl;
    }
    
    if (VERB > 2) {
      printAllWeightsColumns(cerr);
    }
    if (VERB > 3) {
      printAllRawWeightsColumns(cerr, pNorm);
    }
    if (foundPositives > 0 && foundPositivesOldOld > 0 && quickValidation_ &&
           (static_cast<double>(foundPositives - foundPositivesOldOld) <= 
           foundPositivesOldOld * requiredIncreaseOver2Iterations_)) {
      if (VERB > 0) {
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
    printAllWeightsColumns(cerr);
  }
  if (VERB == 3) {
    printAllRawWeightsColumns(cerr, pNorm);
  }
  foundPositives = 0;
  for (size_t set = 0; set < numFolds_; ++set) {
    foundPositives += testScores_[set].calcScores(w_[set], testFdr_);
  }
  if (VERB > 0) {
    std::cerr << "Found " << foundPositives << 
                 " test set PSMs with q<" << testFdr_ << "." << std::endl;
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
int CrossValidation::doStep(bool updateDOC, Normalizer* pNorm, double selectionFdr) {
  // Setup
  struct options* pOptions = new options;
  pOptions->lambda = 1.0;
  pOptions->lambda_u = 1.0;
  pOptions->epsilon = EPSILON;
  pOptions->cgitermax = CGITERMAX;
  pOptions->mfnitermax = MFNITERMAX;
  int estTruePos = 0;
  
  // for determining an appropriate positive training set, the decoys+1 in the 
  // FDR estimates is too restrictive for small datasets
  bool skipDecoysPlusOne = true; 
  for (int set = 0; set < numFolds_; ++set) {
    trainScores_[set].calcScores(w_[set], selectionFdr, skipDecoysPlusOne);
  }
  
  if (DataSet::getCalcDoc() && updateDOC) {
  #pragma omp parallel for schedule(dynamic, 1)
    for (int set = 0; set < numFolds_; ++set) {
      trainScores_[set].recalculateDescriptionOfCorrect(selectionFdr);
      //trainScores_[set].setDOCFeatures(pNorm); // this overwrites features of overlapping training folds...
      testScores_[set].getDOC().copyDOCparameters(trainScores_[set].getDOC());
      testScores_[set].setDOCFeatures(pNorm);
    }
  }
  
#pragma omp parallel for schedule(dynamic, 1) ordered
  for (int set = 0; set < numFolds_; ++set) {
    struct vector_double* pWeights = new vector_double;
    pWeights->d = FeatureNames::getNumFeatures() + 1;
    pWeights->vec = new double[pWeights->d];
      
    double bestCpos = 1, bestCfrac = 1;
      
    int estTruePosFold = processSingleFold(set, selectionFdr, candidatesCpos_, 
					   candidatesCfrac_, bestCpos, bestCfrac, 
					   pWeights, pOptions);
#pragma omp critical (add_tps)
    {
      estTruePos += estTruePosFold;
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
 * @param cposCandidates candidate soft margin parameters for positives
 * @param cfracCandidates candidate soft margin parameters for fraction neg/pos
 * @param bestCpos best soft margin parameter for positives
 * @param bestCfrac best soft margin parameter for fraction neg/pos
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
*/
int CrossValidation::processSingleFold(unsigned int set, double selectionFdr,
    const vector<double>& cposCandidates, const vector<double>& cfracCandidates, 
    double &bestCpos, double &bestCfrac, vector_double* pWeights, 
    options * pOptions) {
  // for determining the number of positives, the decoys+1 in the FDR estimates 
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true;
  
  int bestTruePos = 0;
  vector<double> ww = w_[set]; // SVM weights initial guess and result holder
  vector<double> bestW = w_[set]; // SVM weights with highest true pos estimate
  
  if (VERB > 3) {
    cerr << "Starting processing CV split " << set + 1 << " out of "
         << numFolds_ << endl;
  }
  
  std::vector<Scores> nestedTrainScores(nestedXvalBins_, usePi0_), nestedTestScores(nestedXvalBins_, usePi0_);
  if (nestedXvalBins_ > 1) {
    FeatureMemoryPool featurePool;
  #pragma omp ordered
    {
      trainScores_[set].createXvalSetsBySpectrum(nestedTrainScores, nestedTestScores, nestedXvalBins_, featurePool);
    }
  } else {
    // sub-optimal cross validation
    nestedTrainScores[0] = trainScores_[set];
    nestedTestScores[0] = trainScores_[set];
  }
  
  AlgIn* svmInput = svmInputs_[set % numAlgInObjects_];
  unsigned int nestedFold = 0;
  double cpos = 0;
  double cfrac = 0;

  std::map<std::pair<double, double>, int> intermediateResults;
  std::vector<cpCnTriple>::const_iterator itClassweights = classWeightsPerFold_.begin();
  for ( ; itClassweights != classWeightsPerFold_.end(); ++itClassweights) {
    cpCnTriple cpCnFold = *itClassweights;
    if(cpCnFold.set != set){
      continue;
    }
    nestedFold = cpCnFold.nestedSet;
    cpos = cpCnFold.cpos;
    cfrac = cpCnFold.cfrac;
    nestedTrainScores[nestedFold].generateNegativeTrainingSet(*svmInput, 1.0);
    nestedTrainScores[nestedFold].generatePositiveTrainingSet(*svmInput, selectionFdr, 1.0, trainBestPositive_);
    
    if (VERB > 2) {
      cerr << "Split " << set + 1 << ": Training with " 
        << svmInput->positives << " positives and "
        << svmInput->negatives << " negatives" << std::endl;
    }
    
    // Create storage vector for SVM algorithm
    struct vector_double* Outputs = new vector_double;
    size_t numInputs = svmInput->positives + svmInput->negatives;
    Outputs->vec = new double[numInputs];
    Outputs->d = numInputs;

    if (VERB > 3) cerr << "- cross-validation with Cpos=" << cpos
		       << ", Cneg=" << cfrac * cpos << endl;
    int tp = 0;
    for (int ix = 0; ix < pWeights->d; ix++) {
      pWeights->vec[ix] = 0;
    }
    for (int ix = 0; ix < Outputs->d; ix++) {
      Outputs->vec[ix] = 0;
    }
    svmInput->setCost(cpos, cfrac * cpos);
        
    // Call SVM algorithm (see ssl.cpp)
    L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs);
        
    for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
      ww[i] = pWeights->vec[i];
      cpCnFold.ww[i] = pWeights->vec[i];
    }
        
    tp = nestedTestScores[nestedFold].calcScores(ww, testFdr_, skipDecoysPlusOne);
    if (VERB > 3) {
      cerr << "- cross-validation found " << tp
	   << " training set PSMs with q_liberal<" << testFdr_ << "." << endl;
    }
        
    if (nestedXvalBins_ > 1) {
      intermediateResults[std::make_pair(cpos, cfrac)] += tp;
    } else {
      if (tp >= bestTruePos) {
	if (VERB > 3) {
	  cerr << "Better than previous result, store this." << endl;
	}
	bestTruePos = tp;
	bestW = ww;
	bestCpos = cpos;
	bestCfrac = cfrac;
      }
    }
    delete[] Outputs->vec;
    delete Outputs;
  }
  
  if (nestedXvalBins_ > 1) {
    // Find soft margin parameters with highest estimate of true positives
    std::vector<double>::const_iterator itCpos = cposCandidates.begin();
    for ( ; itCpos != cposCandidates.end(); ++itCpos) {
      double cpos = *itCpos;  
      std::vector<double>::const_iterator itCfrac = cfracCandidates.begin();
      for ( ; itCfrac != cfracCandidates.end(); ++itCfrac) {
        double cfrac = *itCfrac;
        int tp = intermediateResults[std::make_pair(cpos, cfrac)];
        if (tp >= bestTruePos) {
          if (VERB > 3) {
            cerr << "Better than previous result, store this: tp = " << tp << ", cpos = " << cpos << ", cneg = " << cfrac*cpos << endl;
          }
          bestTruePos = tp;
          bestCpos = cpos;
          bestCfrac = cfrac;
        }
      }
    }    
    trainScores_[set].generateNegativeTrainingSet(*svmInput, 1.0);
    trainScores_[set].generatePositiveTrainingSet(*svmInput, selectionFdr, 1.0, trainBestPositive_);
    
    // Create storage vector for SVM algorithm
    struct vector_double* Outputs = new vector_double;
    size_t numInputs = svmInput->positives + svmInput->negatives;
    Outputs->vec = new double[numInputs];
    Outputs->d = numInputs;
    
    for (int ix = 0; ix < pWeights->d; ix++) {
      pWeights->vec[ix] = 0;
    }
    for (int ix = 0; ix < Outputs->d; ix++) {
      Outputs->vec[ix] = 0;
    }
    svmInput->setCost(bestCpos, bestCpos * bestCfrac);
    
    // Call SVM algorithm (see ssl.cpp)
    L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs);
    
    for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
      bestW[i] = pWeights->vec[i];
    }
    delete[] Outputs->vec;
    delete Outputs;
  }
  
  bestTruePos = trainScores_[set].calcScores(bestW, testFdr_);
  
  if (VERB > 2) {
    std::cerr << "Split " << set + 1 << ": Found " << 
        bestTruePos << " training set PSMs with q<" << testFdr_ <<
        " for hyperparameters Cpos=" << bestCpos << 
        ", Cneg=" << bestCfrac * bestCpos << "." << std::endl;
  }
  
  w_[set] = bestW;
  return bestTruePos;
}

void CrossValidation::postIterationProcessing(Scores& fullset,
                                              SanityCheck* pCheck) {
  if (!pCheck->validateDirection(w_)) {
    for (int set = 0; set < numFolds_; ++set) {
      testScores_[set].calcScores(w_[0], selectionFdr_);
    }
  }
  if (DataSet::getCalcDoc()) {
    // TODO: take the average instead of the first DOC model?
    fullset.getDOC().copyDOCparameters(testScores_[0].getDOC());
  }
  fullset.merge(testScores_, selectionFdr_);
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printSetWeights(ostream & weightStream, unsigned int set) {
  weightStream << DataSet::getFeatureNames().getFeatureNames() << 
                  "\tm0" << std::endl;
  weightStream.precision(3);
  weightStream << w_[set][0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << w_[set][ix];
  }
  weightStream << endl;
}

/**
 * Prints the weights of the raw vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printRawSetWeights(ostream & weightStream, unsigned int set, 
                                      Normalizer * pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(w_[set], ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

// used to save weights for reuse as weight input file
void CrossValidation::printAllWeights(ostream & weightStream, 
                                      Normalizer * pNorm) {
  for (unsigned int ix = 0; ix < numFolds_; ++ix) {
    printSetWeights(weightStream, ix);
    printRawSetWeights(weightStream, ix, pNorm);
  }
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

void CrossValidation::printAllWeightsColumns(
    std::vector< std::vector<double> > weightMatrix, ostream & outputStream) {
  // write to intermediate stream to prevent the fixed precision from sticking
  ostringstream weightStream;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    weightStream << " Split" << set + 1 << '\t'; // right-align with weights
  }
  weightStream << "FeatureName" << std::endl;
  size_t numRows = FeatureNames::getNumFeatures() + 1;
  for (unsigned int ix = 0; ix < numRows; ix++) {
    for (unsigned int set = 0; set < numFolds_; ++set) {
      // align positive and negative weights
      if (weightMatrix[set][ix] >= 0) weightStream << " ";
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

void CrossValidation::printAllWeightsColumns(ostream & weightStream) {
  std::cerr << "Learned normalized SVM weights for the " << numFolds_ 
            << " cross-validation splits:" << std::endl;
  printAllWeightsColumns(w_, weightStream);
}

void CrossValidation::printAllRawWeightsColumns(ostream & weightStream, 
                                      Normalizer * pNorm) {
  std::cerr << "Learned raw SVM weights for the " << numFolds_ 
            << " cross-validation splits:" << std::endl;
  std::vector< std::vector<double> > ww = w_;
  for (size_t set = 0; set < w_.size(); ++set) {
    pNorm->unnormalizeweight(w_[set], ww[set]);
  }
  printAllWeightsColumns(ww, weightStream);
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
