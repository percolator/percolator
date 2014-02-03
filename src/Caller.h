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
#include <ctime>
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
#include <boost/foreach.hpp>
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif
#include "Globals.h"
#include "Option.h"
#include "SanityCheck.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "SetHandler.h"
#include "XMLInterface.h"
#include "ssl.h"
#include "MassHandler.h"
#include "ProteinProbEstimator.h"
#include "FidoInterface.h"
#include "CrossValidation.h"

class Caller {
  public:
    enum SetHandlerType {
      NORMAL = 1, SHUFFLED = -1, SHUFFLED_TEST = 2, SHUFFLED_THRESHOLD = 3
    };
    
  public:
    
    Caller();
    Caller(bool uniquePeptides);
    virtual ~Caller();
    
    static string greeter();
    string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    int readFiles();
    
    void fillFeatureSets();
    void calculatePSMProb(bool uniquePeptideRun, Scores *fullset, time_t& procStart,
        clock_t& procStartClock, double& diff, bool TDC = false);
    
    void calculateProteinProbabilitiesFido();
    
    int run();
    
  protected:
    
    XMLInterface xmlInterface;
    
    Normalizer * pNorm;
    SanityCheck * pCheck;
    ProteinProbEstimator* protEstimator;
    
    CrossValidation crossValidation;
    
    string forwardTabInputFN;
    string resultFN;
    string tabFN;
    
    char* xmlInputDir;
    string weightFN;
    
    bool tabInput;
    bool readStdIn;
    
    bool reportUniquePeptides;
    bool target_decoy_competition;
    
    double test_fdr;
    
    double threshTestRatio;
    double trainRatio;
    
    string call;
    string decoyOut;
    
    time_t startTime;
    clock_t startClock;
    //SetHandler normal, shuffled; //,shuffledTest,shuffledThreshold;
    SetHandler setHandler;
    Scores fullset; //,thresholdset;
    map<int, double> scan2rt;

    
};

#endif /*CALLER_H_*/
