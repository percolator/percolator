/*******************************************************************************
 * Percolator v 1.00
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: IntraSetRelation.cpp,v 1.5 2007/02/04 04:50:39 lukall Exp $
 *******************************************************************************/
#include <set>
#include <map>
#include <string>
using namespace std;
#include "IntraSetRelation.h"

IntraSetRelation::IntraSetRelation()
{
}

IntraSetRelation::~IntraSetRelation()
{
}
void IntraSetRelation::registerRel(string pep, set<string> &prot){
  if(numPeptides.count(pep)==0) numPeptides[pep]=1;
  else numPeptides[pep]++;
  set<string>::iterator iProt;
  for(iProt=prot.begin();iProt!=prot.end();iProt++) {
    if(numProteins.count(*iProt)==0) numProteins[*iProt]=1;
    else numProteins[*iProt]++;
    prot2pep[*iProt].insert(pep);
  }
}

int IntraSetRelation::getPepSites(set<string> &prot) {
  int maxn=1;
  set<string>::iterator iProt;
  for(iProt=prot.begin();iProt!=prot.end();iProt++) {
    if (prot2pep.count(*iProt)!=0) {
  	  int n=prot2pep[*iProt].size();
  	  if (n>maxn) maxn=n;
    }
  }
  return maxn;
}

int IntraSetRelation::getNumProt(set<string> &prot) {
  int maxn=1;
  set<string>::iterator iProt;
  for(iProt=prot.begin();iProt!=prot.end();iProt++) {
    string pr = *iProt;
    int n=getNumProt(pr);
    if (n>maxn) maxn=n;
  }
  return maxn;
}
