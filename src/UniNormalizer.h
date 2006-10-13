#ifndef UNINORMALIZER_H_
#define UNINORMALIZER_H_

class UniNormalizer : public Normalizer // virtual Normalizer
{
public:
	UniNormalizer();
	virtual ~UniNormalizer();
    void normalize(const double *in,double* out);
    void unnormalizeweight(const double *in,double* out);
    void normalizeweight(const double *in,double* out);
    void setSet(SetHandler *);
protected:
	vector<double> sub;
	vector<double> div;
};

#endif /*UNINORMALIZER_H_*/
