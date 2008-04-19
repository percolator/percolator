/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SqtSanityCheck.h,v 1.2 2008/04/19 21:40:33 lukall Exp $
 *******************************************************************************/
#ifndef SQTSANITYCHECK_H_
#define SQTSANITYCHECK_H_

#include "SanityCheck.h"

class SqtSanityCheck : public SanityCheck
{
public:
  SqtSanityCheck();
  virtual ~SqtSanityCheck();
  virtual bool validateDirection(double* w);
protected:
  virtual void getDefaultDirection(double *w);
};

#endif /*SQTSANITYCHECK_H_*/
