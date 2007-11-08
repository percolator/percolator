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
