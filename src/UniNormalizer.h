/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: UniNormalizer.h,v 1.8 2007/02/13 18:17:15 lukall Exp $
 *******************************************************************************/
#ifndef UNINORMALIZER_H_
#define UNINORMALIZER_H_

class UniNormalizer : public Normalizer // virtual Normalizer
{
public:
	UniNormalizer();
	virtual ~UniNormalizer();
    virtual void setSet(set<DataSet *> & setVec);
    void normalize(const double *in,double* out);
    void unnormalizeweight(const double *in,double* out);
    void normalizeweight(const double *in,double* out);
protected:
	vector<double> sub;
	vector<double> div;
};

#endif /*UNINORMALIZER_H_*/
