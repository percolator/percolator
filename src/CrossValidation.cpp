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

const unsigned int CrossValidation::xval_fold = 3; /* number of folds for cross validation*/
const double CrossValidation::requiredIncreaseOver2Iterations = 0.01; /* checks cross validation convergence */

CrossValidation::CrossValidation() :
      quickValidation(false), reportPerformanceEachIteration(false),
      testFdr(0.01), selectionFdr(0.01), selectedCpos(0), selectedCneg(0), niter(10) {}


CrossValidation::~CrossValidation() { 
  if (svmInput) {
    delete svmInput;
  }
  svmInput = NULL;
}

void CrossValidation::printParameters(ostringstream & oss) {
  oss << "Hyperparameters fdr=" << selectionFdr
      << ", Cpos=" << selectedCpos << ", Cneg=" << selectedCneg
      << ", maxNiter=" << niter << endl;
}


/** 
 * Sets up the SVM classifier: 
 * - divide dataset into training and test sets for each fold
 * - set parameters (fdr, soft margin)
 * @param w list of normal vectors
 * @return number of positives for initial setup
 */
int CrossValidation::preIterationSetup(Scores & fullset, SanityCheck * pCheck, Normalizer * pNorm) {
  
  w = vector<vector<double> > (xval_fold,vector<double> (FeatureNames::getNumFeatures()+ 1));
  
  svmInput = new AlgIn(fullset.size(), FeatureNames::getNumFeatures() + 1); // One input set, to be reused multiple times
  assert( svmInput );

  if (selectedCpos >= 0 && selectedCneg >= 0) {
    xv_train.resize(xval_fold);
    xv_test.resize(xval_fold);
    
  	fullset.createXvalSetsBySpectrum(xv_train, xv_test, xval_fold);

    if (selectionFdr <= 0.0) {
      selectionFdr = testFdr;
    }
    if (selectedCpos > 0) {
      xv_cposs.push_back(selectedCpos);
    } else {
      xv_cposs.push_back(10);
      xv_cposs.push_back(1);
      xv_cposs.push_back(0.1);
      if (VERB > 0) {
        cerr << "selecting cpos by cross validation" << endl;
      }
    }
    if (selectedCpos > 0 && selectedCneg > 0) {
      xv_cfracs.push_back(selectedCneg / selectedCpos);
    } else {
      xv_cfracs.push_back(1.0 * fullset.getTargetDecoySizeRatio());
      xv_cfracs.push_back(3.0 * fullset.getTargetDecoySizeRatio());
      xv_cfracs.push_back(10.0 * fullset.getTargetDecoySizeRatio());
      if (VERB > 0) {
        cerr << "selecting cneg by cross validation" << endl;
      }
    }
    return pCheck->getInitDirection(xv_test, xv_train, pNorm, w, testFdr);
  } else {
    vector<Scores> myset(1, fullset);
    cerr << "B" << endl;
    return pCheck->getInitDirection(myset, myset, pNorm, w, testFdr);
  }
}

/** 
 * Train the SVM using several cross validation iterations
 * @param w list of normal vectors
 */
void CrossValidation::train(Normalizer * pNorm) {

  if (VERB > 0) {
    cerr << "---Training with Cpos";
    if (selectedCpos > 0) {
      cerr << "=" << selectedCpos;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", Cneg";
    if (selectedCneg > 0) {
      cerr << "=" << selectedCneg;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", fdr=" << selectionFdr << endl;
  }
  
  // iterate
  int foundPositivesOldOld=0, foundPositivesOld=0, foundPositives=0; 
  for (unsigned int i = 0; i < niter; i++) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << " :\t";
    }
    
    foundPositives = xv_step(true);
    
    if (VERB > 1) {
      cerr << "After the iteration step, " << foundPositives
          << " target PSMs with q<" << selectionFdr
          << " were estimated by cross validation" << endl;
    }
    if (VERB > 2) {
      cerr << "Obtained weights" << endl;
      printAllWeights(cerr, pNorm);
    }
    if (foundPositives>0 && foundPositivesOldOld>0 && quickValidation) {
      if ((double)(foundPositives-foundPositivesOldOld)<=(foundPositivesOldOld*requiredIncreaseOver2Iterations)) {
        if (VERB > 1) {
          cerr << "Performance increase over the last two iterations indicate that the algorithm has converged" << endl;
          cerr << "(" << foundPositives << " vs " << foundPositivesOldOld << ")" << endl;
        }
        break;
      }
    }    
    foundPositivesOldOld=foundPositivesOld;    
    foundPositivesOld=foundPositives;
  }
  if (VERB == 2) {
    cerr
    << "Obtained weights (only showing weights of first cross validation set)"
    << endl;
    printSetWeights(cerr, 0, pNorm);
  }
  foundPositives = 0;
  for (size_t set = 0; set < xval_fold; ++set) {
    if (DataSet::getCalcDoc()) {
      xv_test[set].getDOC().copyDOCparameters(xv_train[set].getDOC());
      xv_test[set].setDOCFeatures();
    }
    foundPositives += xv_test[set].calcScores(w[set], testFdr);
  }
  if (VERB > 0) {
    cerr << "After all training done, " << foundPositives << " target PSMs with q<"
        << testFdr << " were found when measuring on the test set"
        << endl;
  }  
}


/** 
 * Executes a cross validation step
 * @param w list of the bins' normal vectors (in linear algebra sense) of the hyperplane from SVM
 * @param updateDOC boolean deciding to recalculate retention features @see DescriptionOfCorrect
 * @return Estimation of number of true positives
 */
int CrossValidation::xv_step(bool updateDOC) {
  // Setup
  struct options* pOptions = new options;
  pOptions->lambda = 1.0;
  pOptions->lambda_u = 1.0;
  pOptions->epsilon = EPSILON;
  pOptions->cgitermax = CGITERMAX;
  pOptions->mfnitermax = MFNITERMAX;
  struct vector_double* pWeights = new vector_double;
  pWeights->d = FeatureNames::getNumFeatures() + 1;
  pWeights->vec = new double[pWeights->d];
  int estTP = 0;
  double best_cpos = 1, best_cfrac = 1;
  if (!quickValidation) {
    for (unsigned int set = 0; set < xval_fold; ++set) {
      estTP += xv_process_one_bin(set,updateDOC, xv_cposs, xv_cfracs, best_cpos, best_cfrac, pWeights, pOptions);   
    }
  } else {
    // Use limited internal cross validation, i.e take the cpos and cfrac values of the first bin 
    // and use it for the subsequent bins 
    estTP += xv_process_one_bin(0,updateDOC, xv_cposs, xv_cfracs, best_cpos, best_cfrac, pWeights, pOptions);
    vector<double> cp(1),cf(1);
    cp[0]=best_cpos; cf[0]= best_cfrac;
    for (unsigned int set = 1; set < xval_fold; ++set) {
      estTP += xv_process_one_bin(set,updateDOC, cp, cf, best_cpos, best_cfrac, pWeights, pOptions);   
    }
  }
  delete[] pWeights->vec;
  delete pWeights;
  delete pOptions;
  return estTP / (xval_fold - 1);
}

/** 
 * Train one of the crossvalidation bins 
 * @param set identification number of the bin that is processed
 * @param w list of normal vectors (in the linear algebra sense) of the hyperplane from SVM, one for each bin
 * @param updateDOC boolean deciding to calculate retention features @see DescriptionOfCorrect
 * @param cpos_vec vector with soft margin parameter for positives
 * @param cfrac_vec vector with soft margin parameter for fraction negatives / positives
 * @param best_cpos best soft margin parameter for positives
 * @param best_cfrac best soft margin parameter for fraction negatives / positives
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
*/
int CrossValidation::xv_process_one_bin(unsigned int set, bool updateDOC, vector<double>& cpos_vec, 
                               vector<double>& cfrac_vec, double &best_cpos, double &best_cfrac, vector_double* pWeights,
                               options * pOptions) {
  int bestTP = 0;
  if (VERB > 2) {
    cerr << "cross validation - fold " << set + 1 << " out of "
         << xval_fold << endl;
  }
  
  vector<double> ww = w[set]; // normal vector initial guess and result holder
  vector<double> bestW = w[set]; // normal vector with highest true positive estimate
  xv_train[set].calcScores(ww, selectionFdr);
  if (DataSet::getCalcDoc() && updateDOC) {
    xv_train[set].recalculateDescriptionOfCorrect(selectionFdr);
  }
  xv_train[set].generateNegativeTrainingSet(*svmInput, 1.0);
  xv_train[set].generatePositiveTrainingSet(*svmInput, selectionFdr, 1.0);
  if (VERB > 2) {
    cerr << "Calling with " << svmInput->positives << " positives and "
         << svmInput->negatives << " negatives\n";
  }
  
  // Create storage vector for SVM algorithm
  struct vector_double* Outputs = new vector_double;
  Outputs->vec = new double[svmInput->positives + svmInput->negatives];
  Outputs->d = svmInput->positives + svmInput->negatives;
  
  // Find combination of soft margin parameters with highest estimate of true positives
  BOOST_FOREACH (const double cpos, cpos_vec) {
    BOOST_FOREACH (const double cfrac, cfrac_vec) {
      if (VERB > 2) cerr << "-cross validation with cpos=" << cpos
          << ", cfrac=" << cfrac << endl;
      int tp = 0;
      for (int ix = 0; ix < pWeights->d; ix++) {
        pWeights->vec[ix] = 0;
      }
      for (int ix = 0; ix < Outputs->d; ix++) {
        Outputs->vec[ix] = 0;
      }
      svmInput->setCost(cpos, (cpos) * (cfrac));
      
      // Call SVM algorithm (see ssl.cpp)
      L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs);
      
      for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
        ww[i] = pWeights->vec[i];
      }
      // sub-optimal cross validation (better would be a set disjoint of the training set)
      tp = xv_train[set].calcScores(ww, testFdr);
      if (VERB > 2) {
        cerr << "- cross validation estimates " << tp
             << " target PSMs over " << testFdr * 100 << "% FDR level"
             << endl;
      }
      if (tp >= bestTP) {
        if (VERB > 2) {
          cerr << "Better than previous result, store this" << endl;
        }
        bestTP = tp;
        bestW = ww;
        best_cpos = cpos;
        best_cfrac = cfrac;
      }
    }
    if (VERB > 2) cerr << "cross validation estimates " << bestTP
        / (xval_fold - 1) << " target PSMs with q<" << testFdr
        << " for hyperparameters Cpos=" << best_cpos << ", Cneg="
        << best_cfrac * best_cpos << endl;
  }
  w[set]=bestW;
  delete[] Outputs->vec;
  delete Outputs;
  return bestTP;
}

void CrossValidation::postIterationProcessing(Scores & fullset, SanityCheck * pCheck) {
  if (!pCheck->validateDirection(w)) {
    fullset.calcScores(w[0]);
  }
  if (VERB > 0) {
    cerr << "Merging results from " << xv_test.size() << " datasets" << endl;
  }
  fullset.merge(xv_test, selectionFdr);
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w normal vector
 */
void CrossValidation::printSetWeights(ostream & weightStream, unsigned int set, Normalizer * pNorm) {
  weightStream << "# first line contains normalized weights, second line the raw weights" << endl;
  weightStream << DataSet::getFeatureNames().getFeatureNames() << "\tm0" << endl;
  weightStream.precision(3);
  weightStream << w[set][0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << w[set][ix];
  }
  weightStream << endl;
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(w[set], ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

void CrossValidation::printAllWeights(ostream & weightStream, Normalizer * pNorm) {
  for (unsigned int ix = 0; ix < xval_fold; ++ix) {
    printSetWeights(weightStream, ix, pNorm);
  }
}

void CrossValidation::printDOC() {
  cerr << "For the cross validation sets the average deltaMass are ";
  for (size_t ix = 0; ix < xv_test.size(); ix++) {
    cerr << xv_test[ix].getDOC().getAvgDeltaMass() << " ";
  }
  cerr << "and average pI are ";
  for (size_t ix = 0; ix < xv_test.size(); ix++) {
    cerr << xv_test[ix].getDOC().getAvgPI() << " ";
  }
  cerr << endl;
}
