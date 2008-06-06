/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SqtSanityCheck.h,v 1.4 2008/06/06 17:13:32 lukall Exp $
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
