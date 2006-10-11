#ifndef STDVNORMALIZER_H_
#define STDVNORMALIZER_H_

class StdvNormalizer : public Normalizer // virtual Normalizer
{
public:
	StdvNormalizer();
	virtual ~StdvNormalizer();
    void normalize(const double *in,double* out);
    void unnormalizeweight(const double *in,double* out);
    void normalizeweight(const double *in,double* out);
    void setSet(IsoChargeSet *);
protected:
	vector<double> avg;
	vector<double> stdv;
};

#endif /*STDVNORMALIZER_H_*/
