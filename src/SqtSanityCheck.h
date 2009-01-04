/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: SqtSanityCheck.h,v 1.5 2009/01/04 22:49:30 lukall Exp $
 *******************************************************************************/
#ifndef SQTSANITYCHECK_H_
#define SQTSANITYCHECK_H_

#include "SanityCheck.h"

class SqtSanityCheck : public SanityCheck
{
public:
  SqtSanityCheck();
  virtual ~SqtSanityCheck();
  virtual bool validateDirection(vector<vector<double> >& w);
protected:
  virtual void getDefaultDirection(vector<vector<double> >& w);
};

#endif /*SQTSANITYCHECK_H_*/
