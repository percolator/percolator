/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Normalizer.h,v 1.10 2007/02/13 18:17:15 lukall Exp $
 *******************************************************************************/
#ifndef NORMALIZER_H_
#define NORMALIZER_H_

class Normalizer
{
public:
	virtual ~Normalizer();
    virtual void setSet(set<DataSet *> & setVec){;}
    void normalizeSet(set<DataSet *> & setVec);
    virtual void normalize(const double *in,double* out){;}
    virtual void unnormalizeweight(const double *in,double* out){;}
    virtual void normalizeweight(const double *in,double* out){;}
    static Normalizer * getNew();
    static void setType(int type);
	const static int UNI = 0;
	const static int STDV = 1;
protected:
	Normalizer();
	static int subclass_type;
};

#endif /*NORMALIZER_H_*/
