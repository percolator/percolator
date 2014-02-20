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
#include "ssl.h"


class CrossValidation {
  
  public:
    CrossValidation();
    ~CrossValidation();
    
    int preIterationSetup(Scores & fullset, SanityCheck * pCheck, Normalizer * pNorm);
    
    void train(Normalizer * pNorm);
    
    void postIterationProcessing(Scores & fullset, SanityCheck * pCheck);
    
    void printParameters(ostringstream & oss);
    void printSetWeights(ostream & weightStream, unsigned int set, Normalizer * pNorm);
    void printAllWeights(ostream & weightStream, Normalizer * pNorm);
    void printDOC();
    
    void inline setSelectedCpos(double cpos) { selectedCpos = cpos; }
    double inline getSelectedCpos() { return selectedCpos; }
    void inline setSelectedCneg(double cneg) { selectedCneg = cneg; }
    double inline getSelectedCneg() { return selectedCneg; }
    void inline setSelectionFdr(double fdr) { selectionFdr = fdr; }
    double inline getSelectionFdr() { return selectionFdr; }
    
    void inline setTestFdr(double fdr) { testFdr = fdr; }
    void inline setQuickValidation(bool on) { quickValidation = on; }
    void inline setReportPerformanceEachIteration(bool on) { reportPerformanceEachIteration = on; }
    unsigned int getNiter() { return niter; }
    
  protected:
    AlgIn *svmInput;
    vector< vector<double> > w;
    
    bool quickValidation;
    bool reportPerformanceEachIteration;
    
    double testFdr;
    double selectionFdr;
    double selectedCpos;
    double selectedCneg;
    
    unsigned int niter;
    const static double requiredIncreaseOver2Iterations;
    
    const static unsigned int xval_fold;
    vector<Scores> xv_train, xv_test;
    vector<double> xv_cposs, xv_cfracs;
    
    int xv_process_one_bin(unsigned int set, bool updateDOC, 
                           vector<double>& cpos_vec, vector<double>& cfrac_vec, 
                           double& best_cpos, double& best_cfrac, 
                           vector_double* pWeights, options * pOptions);
    int xv_step(bool updateDOC = false);
    
};

#endif /*CROSSVALIDATION_H_*/
