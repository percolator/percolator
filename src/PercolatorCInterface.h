#ifndef PERCOLATOR_C_INTERFACE_H_
#define PERCOLATOR_C_INTERFACE_H_
#ifdef __cplusplus
extern "C" {
#endif

/** Number of target and decoy sets that external program will hand over to percolator.
 * The value should correspond to the number of sequence databases that have been searched.
 * Percolators validation strategy will be the same as for the stand alone version given the
 * corresponding number of sqt-files as input. */
typedef enum {TWO_SETS=2,THREE_SETS,FOUR_SETS} NSet;

/** Identifying which set the PSM belongs to*/
typedef enum {TARGET=0,DECOY1,DECOY2,DECOY3} SetType;

/** Call that initiates percolator */
void pcInitiate(NSet sets, unsigned int numFeatures, unsigned int numSpectra, char ** featureNames, double pi0);

/** Call that sets verbosity level
 *  0 is quiet, 2 is default, 5 is more than you want */
void pcSetVerbosity(int verbosity);


/** Register a PSM */
void pcRegisterPSM(SetType set, char * identifier, double * features);

/** Function called when we want to start processing */
void pcExecute(); 

/**
 * Given the set enum and features, return the Percolator score for the PSM
 */
void pcScorePSM(
  SetType set, ///< The PSM tag -in
  double* features,  ///< the features -in
  double* score ///< output the Percolator score -out
  );

/** Function called when retrieving target scores and q-values after processing,
  * the array should be numSpectra long and will be filled in the same order
  * as the features were inserted */
void pcGetScores(double *scoreArr, double *qArr); 

/** Function that should be called after processing finished */
void pcCleanUp(); 

#ifdef __cplusplus
}
#endif
#endif /*PERCOLATOR_C_INTERFACE_H_*/
