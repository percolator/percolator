#ifndef NORMALIZER_H_
#define NORMALIZER_H_

#include <vector>
#include "DataSet.h"
#include "IsoChargeSet.h"

class Normalizer
{
public:
	Normalizer();
	virtual ~Normalizer();
    void setSets(IsoChargeSet *set);
    void normalize(double *in,double* out);
protected:
	vector<double> avg;
	vector<double> stdv;
};

#endif /*NORMALIZER_H_*/
