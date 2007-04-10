#ifndef REMOTECONTROL_H_
#define REMOTECONTROL_H_
#ifdef __cplusplus
extern "C" {
#endif

/* Number of target and decoy sets that external program will hand over to percolator.
 * The value should correspond to the number of sequence databases that have been searched.
 * Percolators validation strategy will be the same as for the stand alone version given the
 * corresponding number of sqt-files as input. */
typedef enum {TWO_SETS=2,THREE_SETS,FOUR_SETS} NSet;

/* Identifying which set the PSM belongs to*/
typedef enum {TARGET=0,DECOY1,DECOY2,DECOY3} SetType;

/* Call that initiates percolator */
void initiate(NSet sets, unsigned int numOfFeatures, unsigned int numOfSpectra, char ** featureNames);

/* Register a PSM */
void registerPSM(SetType set, char * identifier, double * features);

/* Function called after all features are registered */
void endRegistration();

/* Function called when we want to start processing */
void go(); 

#ifdef __cplusplus
}
#endif
#endif /*REMOTECONTROL_H_*/
