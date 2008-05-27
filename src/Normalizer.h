/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Normalizer.h,v 1.12 2008/05/27 23:09:08 lukall Exp $
 *******************************************************************************/
#ifndef NORMALIZER_H_
#define NORMALIZER_H_
class DataSet;
#include <set>
using namespace std;

class Normalizer
{
public:
	virtual ~Normalizer();
    virtual void setSet(set<DataSet *> & setVec){;}
    void normalizeSet(set<DataSet *> & setVec);
    virtual void normalize(const double * in, double * out){;}
    virtual void unnormalizeweight(const vector<double>& in,vector<double>& out){;}
    virtual void normalizeweight(const vector<double>& in,vector<double>& out){;}
    static Normalizer * getNew();
    static void setType(int type);
	const static int UNI = 0;
	const static int STDV = 1;
protected:
	Normalizer();
	static int subclass_type;
};

#endif /*NORMALIZER_H_*/
