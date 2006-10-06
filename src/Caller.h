#ifndef CALLER_H_
#define CALLER_H_
#include "Normalizer.h"
#include "IsoChargeSet.h"
#include "Scores.h"
class Caller
{
public:
	Caller();
	virtual ~Caller();
	inline void setSet(IsoChargeSet *p){pSet=p;}
	void step(double *);
protected:
    IsoChargeSet *pSet;
    Scores scores;
};

#endif /*CALLER_H_*/
