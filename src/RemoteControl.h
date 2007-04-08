#ifndef REMOTECONTROL_H_
#define REMOTECONTROL_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {TWO_SETS,THREE_SETS,FOUR_SETS} NSet;
typedef enum {TARGET,DECOY,TEST_DECOY,THRESHOLD_DECOY} SetType;

/* Call that initiates percolator */
void initiate(NSet sets, unsigned int numOfFeatures, unsigned int numOfSpectra, char ** featureNames);

/* Register a PSM */
void registerPSM(SetType set, char * identifier, double * features);

/* Function called after all features are registered */
void endRegistration();

/* Function called when we whant to start processing */
void go(); 

#ifdef __cplusplus
}
#endif



#endif /*REMOTECONTROL_H_*/
