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
#include "SanityCheck.h"
using namespace std;

class Caller {
  public:
    enum SetHandlerType {
      NORMAL = 0, SHUFFLED, SHUFFLED_TEST, SHUFFLED_THRESHOLD
    };
  public:
    Caller();
    virtual ~Caller();
    void readRetentionTime(string filename);
    void train(vector<vector<double> >& w);
    int xv_step(vector<vector<double> >& w, bool updateDOC = false);
    static string greeter();
    string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, vector<double>& w);
    void readWeights(istream & weightStream, vector<double>& w);
    void readFiles();
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
    void writeXML(ostream& os, Scores& fullset);
  protected:
    void countTargetsAndDecoys( std::string & fname, unsigned int & nrTargets , unsigned int & nrDecoys );
    Normalizer * pNorm;
    SanityCheck * pCheck;
    AlgIn *svmInput;
    string tokyoCabinetTmpFN;
    string xmlInputFN;
    string xmlOutputFN;
    string modifiedFN;
    string modifiedDecoyFN;
    string forwardFN;
    string decoyFN;
    string decoyWC;
    string resultFN;
    string gistFN;
    string tabFN;
    string xmloutFN;
    string weightFN;
    string call;
    string spectrumFile;
    string decoyOut;
    bool gistInput;
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
    Scores fullset; //,thresholdset;
    map<int, double> scan2rt;
};

#endif /*CALLER_H_*/
