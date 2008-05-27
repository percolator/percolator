/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SqtSanityCheck.h,v 1.3 2008/05/27 23:09:08 lukall Exp $
 *******************************************************************************/
#ifndef SQTSANITYCHECK_H_
#define SQTSANITYCHECK_H_

#include "SanityCheck.h"

class SqtSanityCheck : public SanityCheck
{
public:
  SqtSanityCheck();
  virtual ~SqtSanityCheck();
  virtual bool validateDirection(vector<double>& w);
protected:
  virtual void getDefaultDirection(vector<double>& w);
};

#endif /*SQTSANITYCHECK_H_*/
