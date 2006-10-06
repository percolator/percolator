#ifndef NORMALIZER_H_
#define NORMALIZER_H_

#include <vector>

class IsoChargeSet;

class Normalizer
{
public:
	Normalizer();
	virtual ~Normalizer();
    void setSet(IsoChargeSet *);
    void normalize(const double *in,double* out);
    void unnormalizeweight(const double *in,double* out);
protected:
	vector<double> avg;
	vector<double> stdv;
};

#endif /*NORMALIZER_H_*/
