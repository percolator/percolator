/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SqtSanityCheck.cpp,v 1.5 2008/06/20 23:55:34 lukall Exp $
 *******************************************************************************/
#include "DataSet.h"
#include "Scores.h"
#include "Globals.h"
#include "SanityCheck.h"
#include "SqtSanityCheck.h"

SqtSanityCheck::SqtSanityCheck()
{
}

SqtSanityCheck::~SqtSanityCheck()
{
}

void SqtSanityCheck::getDefaultDirection(vector<vector<double> >& w) {
  // Set init direction to be the most discriminative direction
  for (size_t set = 0; set < w.size();++set) {  
    for (unsigned int ix=0;ix < FeatureNames::getNumFeatures()+1;++ix)
      w[set][ix]=0;
    w[set][2]=1.61;    // deltCn 
    w[set][3]=1.1;     // Xcorr
    w[set][7]=-0.573;  // Peptide length
    w[set][8]=0.0335;  // Charge 1
    w[set][9]=0.149;   // Charge 2
    w[set][10]=-0.156; // Charge 3
  }
/*

# first line contains normalized weights, second line the raw weights
deltCn  Xcorr   PepLen  Charge1 Charge2 Charge3 m0
1.61    1.1     -0.573  0.0335  0.149   -0.156  -5.23
20.9    1.99    -0.105  0.328   0.299   -0.312  -9.06

        feat[0]=log(rSp);                     // rank by Sp
        feat[1]=0.0;                     // delt5Cn (leave until last M line)
        feat[2]=0.0;                     // deltCn (leave until next M line)
        feat[3]=xcorr;                   // Xcorr
        feat[4]=sp;                      // Sp
        feat[5]=matched/expected;        // Fraction matched/expected ions
        feat[6]=mass;                    // Observed mass
        feat[7]=peptideLength(pep);      // Peptide length
        feat[8]=(charge==1?1.0:0.0);     // Charge
        feat[9]=(charge==2?1.0:0.0);
        feat[10]=(charge==3?1.0:0.0); 
*/
  return;
}


bool SqtSanityCheck::validateDirection(vector<vector<double> >& w) {
  bool ok=SanityCheck::validateDirection(w);
  for (size_t set = 0; set < w.size();++set) {
    if (w[set][3]<0) {
      ok=false;
      if (VERB>1) cerr << "Warning weight for XCorr negative" << endl;
    }
    if (w[set][2]<0) {
      ok=false;
      if (VERB>1) cerr << "Warning weight for deltaCn negative" << endl;
    }
  }
  if (!ok)
    resetDirection(w);
  return ok;
}

