#ifndef RESULTHOLDER_H_
#define RESULTHOLDER_H_

class ResultHolder
{
public:
	ResultHolder();
	ResultHolder(const double score,const double q,const double po,const string& i,const string& pe,const string& p);
	virtual ~ResultHolder();
    double score,q,posterior;
    string id,pepSeq,prot;
};

bool operator>(const ResultHolder &one, const ResultHolder &other);
bool operator<(const ResultHolder &one, const ResultHolder &other); 
ostream& operator<<(ostream &out,const ResultHolder &obj);

#endif /*RESULTHOLDER_H_*/
