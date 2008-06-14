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
  Unspecified
};

enum MSFileFormat {
	bms1,
	bms2,
	cms1,
	cms2,
	ms1,
	ms2,
	mzXML,
  mzData,
	zs,
	uzs,
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
  na
};

struct MSHeader {
	char header[16][128];
};

struct MSScanInfo {
	int scanNumber[2];
	int numDataPoints;
	int numZStates;
  float rTime;
  double mz;
};

struct Peak_T {
  double mz;
  float intensity;
};

struct ZState {
  int z;
  double mz;
};

#endif


