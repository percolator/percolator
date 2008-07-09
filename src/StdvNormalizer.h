/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: StdvNormalizer.h,v 1.10 2008/07/09 00:54:19 lukall Exp $
 *******************************************************************************/
#ifndef STDVNORMALIZER_H_
#define STDVNORMALIZER_H_

class StdvNormalizer : public Normalizer // virtual Normalizer
{
public:
	StdvNormalizer();
	virtual ~StdvNormalizer();
    virtual void setSet(set<DataSet *> & setVec, size_t numFeatures, size_t numRetentionFeatures);
    void unnormalizeweight(const vector<double>& in,vector<double>& out);
    void normalizeweight(const vector<double>& in,vector<double>& out);
};

#endif /*STDVNORMALIZER_H_*/
