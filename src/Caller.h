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

#ifndef CALLER_H_
#define CALLER_H_
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <memory>
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif
#include "Option.h"
#include "SanityCheck.h"
#include "SqtSanityCheck.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "SetHandler.h"
#include "XMLInterface.h"
#include "ssl.h"
#include "Globals.h"
#include "MassHandler.h"
#include "ProteinProbEstimator.h"
#include "FidoInterface.h"

class Caller {
  public:
    enum SetHandlerType {
      NORMAL = 1, SHUFFLED = -1, SHUFFLED_TEST = 2, SHUFFLED_THRESHOLD = 3
    };
    
  public:
    
    Caller();
    Caller(bool uniquePeptides);
    virtual ~Caller();
    void train(vector<vector<double> >& w);
    int xv_process_one_bin(unsigned int set, vector<vector<double> >& w, 
                           bool updateDOC, vector<double>& cpos_vec, 
                           vector<double>& cfrac_vec, double& best_cpos, 
                           double &best_cfrac, vector_double* pWeights,
                           options * pOptions);
    int xv_step(vector<vector<double> >& w, bool updateDOC = false);
    static string greeter();
    string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, vector<double>& w);
    void readWeights(istream & weightStream, vector<double>& w);
    int readFiles();
    
    void fillFeatureSets();
    int preIterationSetup(vector<vector<double> >& w);
    void calculatePSMProb(bool uniquePeptideRun, Scores *fullset, time_t& procStart,
        clock_t& procStartClock, vector<vector<double> >& w, double& diff, bool TDC = false);
    
    void calculateProteinProbabilitiesFido();
    
    int run();
    
  protected:
    
    XMLInterface xmlInterface;
    
    Normalizer * pNorm;
    SanityCheck * pCheck;
    AlgIn *svmInput;
    ProteinProbEstimator* protEstimator;
    
    string forwardTabInputFN;
    string decoyWC;
    string resultFN;
    string tabFN;
    
    char* xmlInputDir;
    string weightFN;
    
    bool tabInput;
    bool readStdIn;
    
    bool quickValidation;
    bool reportPerformanceEachIteration;
    
    bool reportUniquePeptides;
    
    bool target_decoy_competition;
    
    double test_fdr;
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    
    double threshTestRatio;
    double trainRatio;
    unsigned int niter;
    
    string call;
    string otherCall;
    string decoyOut;
    bool outputAll;
    bool showExpMass;
    
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    vector<Scores> xv_train, xv_test;
    vector<double> xv_cposs, xv_cfracs;
    //SetHandler normal, shuffled; //,shuffledTest,shuffledThreshold;
    SetHandler setHandler;
    Scores fullset; //,thresholdset;
    map<int, double> scan2rt;
    /*fido parameters*/
    double fido_alpha;
    double fido_beta;
    double fido_gamma;
    bool fido_nogrouProteins; 
    bool fido_trivialGrouping;
    bool fido_noprune;
    bool fido_noseparate;
    bool fido_reduceTree;
    bool fido_truncate;
    unsigned fido_depth;
    double fido_mse_threshold;
    /* general protein probabilities options */
    bool tiesAsOneProtein;
    bool usePi0;
    bool outputEmpirQVal;
    std::string decoy_prefix;

    
};

#endif /*CALLER_H_*/
