#ifndef _MSTOOLKITTYPES_H
#define _MSTOOLKITTYPES_H

enum MSSpectrumType {
  MS1,
  MS2,
  MS3,
  ZS,
  UZS,
  IonSpec,
  SRM,
  REFERENCE,
  Unspecified
};

enum MSFileFormat {
  bms1,
  bms2,
  cms1,
  cms2,
  mgf,
  ms1,
  ms2,
  msmat_ff,
  mzXML,
  mzData,
  raw,
  sqlite,
  psm,
  uzs,
  zs,
  dunno
};

enum MSTag {
  no,
  D,
  H,
  I,
  S,
  Z
};

enum MSActivation {
  CID,
  ECD,
  ETD,
  PQD,
  HCD,
  na
};

struct MSHeader {
	char header[16][128];
};

struct MSScanInfo {
	int scanNumber[2];
	int numDataPoints;
  int numEZStates;
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

struct DataPeak {
	double dMass;
	double dIntensity;
}; //For RAW files

struct Peak_T {
  double mz;
  float intensity;
};

struct ZState {
  int z;
  double mz;   //M+H, not mz
};

struct EZState {
  int z;
  double mh;      //M+H
  float pRTime;   //precursor area
  float pArea;    //precursor retention time
};

#endif


