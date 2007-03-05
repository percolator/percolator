#ifndef SEQUESTOUT_H_
#define SEQUESTOUT_H_

class SequestOut
{
public:
	SequestOut();
	virtual ~SequestOut();
   char cAA1;
   char cAA2;
   string szFileName;
   string szBaseFileName;
   string szProt;
   string szPlainPep;
   string szSubPep;
   string szDSite;
   string szMod;
   string szDup;
   string szDatabase;
   double dAMass;
   double dMass;
   double dXC;
   double dDeltCn;
   double dSp;
   double dMass1;
   double dMass2;
   double dMass3;
   double dMass4;
   double dMass5;
   double dMass6;
   int  iRankSp;
   int  iMassType;
   int  iIon;
   int  iTot;
   int  bSpecialDeltCn;
   int  bNucDb;
};

class Header
{
public:
	Header() {}
	virtual ~Header() {}
    string szDate;
    string szTime;
    string szTimeSuffix;
    string szMassType;

};

#endif /*SEQUESTOUT_H_*/
