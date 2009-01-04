/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: ResultHolder.h,v 1.4 2009/01/04 22:49:30 lukall Exp $
 *******************************************************************************/
#ifndef RESULTHOLDER_H_
#define RESULTHOLDER_H_

class ResultHolder
{
public:
	ResultHolder();
	ResultHolder(const double score,const double q,const double po,const string& i,const string& pe = "",const string& p  = "");
	virtual ~ResultHolder();
    double score,q,posterior;
    string id,pepSeq,prot;
};

bool operator>(const ResultHolder &one, const ResultHolder &other);
bool operator<(const ResultHolder &one, const ResultHolder &other); 
ostream& operator<<(ostream &out,const ResultHolder &obj);

#endif /*RESULTHOLDER_H_*/
