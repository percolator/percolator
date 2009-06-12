#include <assert.h>
#include <cmath>
#include "MassHandler.h"

MassHandler::MassHandler()
{
}

MassHandler::~MassHandler()
{
}

bool MassHandler::monoisotopic=false;

double MassHandler::massDiff(double observedMass, double calculatedMass, unsigned int charge, const string& peptide) {
  assert(charge>0);
  double dm = observedMass-calculatedMass;
  if (monoisotopic) {
    double isodm = dm-1;
    for(int isotope=0;isotope<5;++isotope) {
  	  if (abs(isodm)>abs(dm+isotope))
        isodm = dm+isotope;
    }
    dm=isodm/calculatedMass;
    return dm;
  }
  return dm/charge;
}
