/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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
#include "ssl.h"
#include "Globals.h"
#include "MassHandler.h"
#include "Enzyme.h"
#include "config.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include "parser.hxx"
#include "serializer.hxx"
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include "percolator_in.hxx"
#include "GroupPowerBigraph.h"

class Caller {
  public:
    enum SetHandlerType {
      NORMAL = 0, SHUFFLED, SHUFFLED_TEST, SHUFFLED_THRESHOLD
    };
  public:
    Caller();
    Caller(bool uniquePeptides);
    virtual ~Caller();
    void train(vector<vector<double> >& w);
    int xv_step(vector<vector<double> >& w, bool updateDOC = false);
    static string greeter();
    string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, vector<double>& w);
    void readWeights(istream & weightStream, vector<double>& w);
    void readFiles();
    void initializeSanityCheck();
    void filelessSetup(const unsigned int numFeatures,
                       const unsigned int numSpectra, char** fetureNames,
                       double pi0);
    void fillFeatureSets();
    int preIterationSetup(vector<vector<double> >& w);
    Scores* getFullSet() {
      return &fullset;
    }
    int run();
    SetHandler* getSetHandler(SetHandlerType sh) {
      switch (sh) {
        case NORMAL:
          return &normal;
        case SHUFFLED:
          return &shuffled;
        case SHUFFLED_TEST:
          return NULL;
        case SHUFFLED_THRESHOLD:
          return NULL;
        default:
          return NULL;
      }
    }
    int proteinLevelProbabilities(Scores& fullset);
    void writeXML(bool uniquePeptides);
    string xmlOutputFN;
    Scores fullset; //,thresholdset;
  protected:
    void countTargetsAndDecoys( std::string & fname, unsigned int & nrTargets , unsigned int & nrDecoys );
    Normalizer * pNorm;
    SanityCheck * pCheck;
    AlgIn *svmInput;
    string xmlInputFN;
    bool readStdIn;
    string forwardFN;
    string decoyFN;
    string decoyWC;
    string resultFN;
    string tabFN;
    string weightFN;
    string call;
    string otherCall;
    string decoyOut;
    bool outputAll;
    bool tabInput;
    bool dtaSelect;
    bool docFeatures;
    bool reportPerformanceEachIteration;
    bool reportUniquePeptides;
    double test_fdr;
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    double threshTestRatio;
    double trainRatio;
    unsigned int niter;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    vector<Scores> xv_train, xv_test;
    vector<double> xv_cposs, xv_cfracs;
    SetHandler normal, shuffled; //,shuffledTest,shuffledThreshold;
    map<int, double> scan2rt;
};

#endif /*CALLER_H_*/
