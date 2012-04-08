// Written by Oliver Serang 2009
// see license for more information

#include "ProteinIdentifier.h"

ProteinIdentifier::ProteinIdentifier()
{
  
  //TODO these guys should be parametizable
  
  ProteinThreshold = 1e-2;
  //ProteinThreshold = 1e-5;
  
  //PeptideThreshold = 1e-9;
  PeptideThreshold = 1e-3;
  //PeptideThreshold = 9e-3;
  
  PsmThreshold = 0.0;
  
  //PeptidePrior = 0.5;
  PeptidePrior = 0.1;
  //PeptidePrior = .07384;
}

