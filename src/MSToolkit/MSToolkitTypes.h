#ifndef _MSTOOLKITTYPES_H
#define _MSTOOLKITTYPES_H

enum MSSpectrumType {
  MS1, MS2, MS3, ZS, UZS, IonSpec, SRM, REFERENCE, Unspecified
};

enum MSFileFormat {
  bms1, bms2, cms1, cms2, ms1, ms2, msmat_ff, mzXML, mzData, sqlite, psm,
  zs, uzs, dunno
};

enum MSTag {
  no, D, H, I, S, Z
};

enum MSActivation {
  CID, ECD, ETD, PQD, HCD, na
};

struct MSHeader {
    char header[16][128];
};

struct MSScanInfo {
    int scanNumber[2];
    int numDataPoints;
    int numZStates;
    float rTime;
    float IIT;
    float BPI;
    double mz;
    double convA;
    double convB;
    double TIC;
    double BPM;
};

struct Peak_T {
    double mz;
    float intensity;
};

struct ZState {
    int z;
    double mz; //M+H, not mz
};

#endif

