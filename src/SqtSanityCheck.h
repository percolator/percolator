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
