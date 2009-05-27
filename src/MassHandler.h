#ifndef MASSHANDLER_H_
#define MASSHANDLER_H_

#include<string> 
using namespace std;

class MassHandler
{
public:
	MassHandler();
	virtual ~MassHandler();
	static void setMonoisotopicMass(bool mi) {monoisotopic=mi;}
	static double massDiff(double observedMass, double calculatedMass, unsigned int charge, const string& peptide);
	static bool monoisotopic;
};

#endif /*MASSHANDLER_H_*/
