// Written by Oliver Serang 2009
// see license for more information

#include "ProteinIdentifier.h"

ProteinIdentifier::ProteinIdentifier()
{
  ProteinThreshold = 1e-5;
  //  ProteinThreshold = 0.0;
  
  //  PeptideThreshold = 1e-9;
  PeptideThreshold = 9e-3;

  //  PeptideProphetPrior = .5;
  PeptideProphetPrior = .07384;
}

