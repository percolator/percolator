#ifndef REMOTECONTROL_H_
#define REMOTECONTROL_H_

typedef enum {NO_ENZYME,TRYPSIN,CHYMOTRYPSIN,ELASTASE} Enzyme;
typedef enum {TWO_SETS,THREE_SETS,FOUR_SETS} NSet;
typedef enum {TARGET,DECOY,TEST_DECOY,THRESHOLD_DECOY} SetType;

struct PSM {
  char * identifier;     /* An unique identifier for this PSM. Should be unique to spectra, charge, rank and set type */
  double XCorr;          /* XCorr score */
  double deltaCn;        /* Fractional difference between current and second best XCorr */
  double deltaLCn;       /* Fractional difference between current and fifth best XCorr */
  double Sp;             /* Sp score */
  unsigned int rankSP;   /* Rank when sorted on Sp score */
  unsigned int charge;   /* The charge state for which we processed the spectrum */
  double observedMass;   /* The mass of the peptide, calculated from the precursor m/z */
  double calculatedMass; /* The mass of the peptide, calculated from the peptide sequence */
  double ionFrac;        /* matched ions/expected ions */
  unsigned int numSP;    /* number of peptides matching the mass criterium */
  char * peptideSeq;     /* The peptide sequence flanked by neighbouring amino acids and dots 
                            i.e. I.AMFAK.E */
  unsigned int numProteinIds; /* Number of proteins where the peptide ocurs */
  char ** proteinIds;     /* The protein ids of the proteins where the peptide occurs */
};


class RemoteControl
{
public:
	RemoteControl();
	virtual ~RemoteControl();
    /* Call that initiates percolator */
    static void initiate(NSet sets, unsigned int numOfSpectra,Enzyme enz);
    /* Register a PSM */
    static void registerPSM(SetType set, PSM psm);
    /* Function called after all features are registered */
    static void endRegistration();
    /* Function called when we whant to start processing */
    static void go(); 
};

#endif /*REMOTECONTROL_H_*/
