/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-9 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: ResultHolder.cpp,v 1.5 2009/01/09 14:41:00 lukall Exp $
 *******************************************************************************/
#include <iostream>
#include <string>
using namespace std;
#include "ResultHolder.h"

ResultHolder::ResultHolder():
score(0.0),q(1.0),posterior(1.0),pepSeq(""),prot("")
{
}

ResultHolder::ResultHolder(const double sc,const double qq,const double po,
                           const string& i,const string& pe,const string& p):
score(sc),q(qq),posterior(po),id(i),pepSeq(pe),prot(p)
{
}

ResultHolder::~ResultHolder()
{
}

bool operator>(const ResultHolder &one, const ResultHolder &other) 
{return (one.score>other.score);}

bool operator<(const ResultHolder &one, const ResultHolder &other) 
{return (one.score<other.score);}

ostream& operator<<(ostream &out,const ResultHolder &obj)
{
    out << obj.id << "\t" << obj.score << "\t" << obj.q << "\t";
    out << obj.posterior << "\t" << obj.pepSeq << obj.prot;
    return out;
}
