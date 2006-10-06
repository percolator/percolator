#ifndef ISOCHARGESET_H_
#define ISOCHARGESET_H_
#include "DataSet.h"
class IsoChargeSet
{
public:
	IsoChargeSet(int charge);
	virtual ~IsoChargeSet();
	void setSet(vector<DataSet> *set);
	const double * const getNext(int& ,int& );
	int const getLabel(int *setPos);
	int const inline getSize() {return n_points;}
protected:
    int charge;
    vector<DataSet> *pSet;
    int n_points;
};

#endif /*ISOCHARGESET_H_*/
