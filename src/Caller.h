#ifndef CALLER_H_
#define CALLER_H_

class Caller
{
public:
	Caller();
	virtual ~Caller();
	inline setSet(IsoChargeSet *p){pSet=p;}
	step(double *);
protected:
    IsoChargeSet *pSet;
    Scores scores;
};

#endif /*CALLER_H_*/
