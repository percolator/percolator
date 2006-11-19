#ifndef STDVNORMALIZER_H_
#define STDVNORMALIZER_H_

class StdvNormalizer : public Normalizer // virtual Normalizer
{
public:
	StdvNormalizer();
	virtual ~StdvNormalizer();
    virtual void setSet(set<DataSet *> & setVec);
    void normalize(const double *in,double* out);
    void unnormalizeweight(const double *in,double* out);
    void normalizeweight(const double *in,double* out);
protected:
	vector<double> avg;
	vector<double> stdv;
};

#endif /*STDVNORMALIZER_H_*/
