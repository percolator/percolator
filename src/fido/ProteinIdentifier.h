// Written by Oliver Serang 2009
// see license for more information

#ifndef _ProteinIdentifier_H
#define _ProteinIdentifier_H

#include "StringTable.h"
#include "Array.h"
#include "Set.h"
#include "Matrix.h"
#include "Random.h"
#include "Model.h"
#include "Scores.h"

using namespace std;

class ProteinIdentifier
{
public:
  ProteinIdentifier();
  virtual ~ProteinIdentifier() {}

  friend Scores* operator >>(Scores* fullset, ProteinIdentifier & pi)
  {
    pi.read(fullset);
    return fullset;
  }

  virtual void printProteinWeights() const = 0;

  class FormatException {};
  
  double ProteinThreshold, PeptideThreshold, PsmThreshold;
  
protected:

  virtual void read(Scores* fullset) = 0;
  
  double PeptidePrior;
};

#endif