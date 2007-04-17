#include <iostream>
using namespace std;
#include "PercolatorCInterface.h"


/** Call that initiates percolator */
void pcInitiate(NSet sets, unsigned int numFeatures, unsigned int numSpectra, char ** featureNames, double pi0) {
	cout << "A C++ call" << endl;
}

/** Register a PSM */
void pcRegisterPSM(SetType set, char * identifier, double * features) {;}

/** Function called after all features are registered */
void pcEndRegistration() {;}

/** Function called when we want to start processing */
void pcExecute() {;} 

/** Function called when retrieving target q-values after processing,
  * the array should be numSpectra long and will be filled in the same order
  * as the features were inserted */
void pcGetQ(double *qArr) {;} 

/** Function called when retrieving target scores after processing,
  * the array should be numSpectra long and will be filled in the same order
  * as the features were inserted */
void pcGetScores(double *scoreArr) {;} 

/** Function that should be called after processing finished */
void pcpcCleanUp() {;}
