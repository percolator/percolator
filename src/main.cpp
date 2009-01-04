/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: main.cpp,v 1.3 2009/01/04 22:49:31 lukall Exp $
 *******************************************************************************/
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"

int main(int argc, char **argv){
  Caller *pCaller = new Caller();
  int retVal = -1;
  if(pCaller->parseOptions(argc,argv))
  {
    retVal = pCaller->run();
  }
  delete pCaller;
  Globals::clean();
  return retVal;
}   
