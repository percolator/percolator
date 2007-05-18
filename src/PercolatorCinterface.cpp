#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"
#include "PercolatorCInterface.h"

static Caller *pCaller=NULL;
static NSet nset=FOUR_SETS;

Caller * getCaller() {
    if (pCaller==NULL) {
      cerr << "Object pCaller not properly assigned" << endl;
      exit(-1);
    }
    return pCaller;
} 


/** Call that initiates percolator */
void pcInitiate(NSet sets, unsigned int numFeatures, unsigned int numSpectra, char ** featureNames, double pi0) {
    pCaller=new Caller();
    nset=sets;
    pCaller->filelessSetup((unsigned int) sets, numFeatures, numSpectra);
}

/** Register a PSM */
void pcRegisterPSM(SetType set, char * identifier, double * features) {;}

/** Function called after all features are registered */
void pcEndRegistration() {;}

/** Function called when we want to start processing */
void pcExecute() {
  bool separateShuffledTestSetHandler = nset>TWO_SETS;
  bool separateShuffledThresholdSetHandler = nset==FOUR_SETS;
  pCaller->fillFeatureSets(separateShuffledTestSetHandler,separateShuffledThresholdSetHandler);
  pCaller->preIterationSetup();
} 

/** Function called when retrieving target q-values after processing,
  * the array should be numSpectra long and will be filled in the same order
  * as the features were inserted */
void pcGetQ(double *qArr) {;} 

/** Function called when retrieving target scores after processing,
  * the array should be numSpectra long and will be filled in the same order
  * as the features were inserted */
void pcGetScores(double *scoreArr) {;} 

/** Function that should be called after processing finished */
void pcCleanUp() {
    if (pCaller) {
      delete pCaller;
      pCaller=NULL;
    }
}
