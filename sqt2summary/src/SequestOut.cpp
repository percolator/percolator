#include <string>
using namespace std;
#include "SequestOut.h"

SequestOut::SequestOut() :
 cAA1('\0'),
 cAA2('\0'),
 szFileName(""),
 szBaseFileName(""),
 szProt(""),
 szPlainPep(""),
 szSubPep(""),
 szDSite(""),
 szMod(""),
 szDup(""),
 szDatabase(""),
 dAMass(0.0),
 dMass(0.0),
 dXC(0.0),
 dDeltCn(0.0),
 dSp(0.0),
 dMass1(0.0),
 dMass2(0.0),
 dMass3(0.0),
 dMass4(0.0),
 dMass5(0.0),
 dMass6(0.0),
 iRankSp(0),
 iMassType(0),
 iIon(0),
 iTot(0),
 bSpecialDeltCn(0),
 bNucDb(0)
{
}

SequestOut::~SequestOut()
{
}
