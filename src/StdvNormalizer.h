/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: StdvNormalizer.h,v 1.8 2007/02/13 18:17:14 lukall Exp $
 *******************************************************************************/
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
