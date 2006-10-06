#ifndef ISOCHARGESET_H_
#define ISOCHARGESET_H_
#include "DataSet.h"
#include "Normalizer.h"
class IsoChargeSet
{
public:
	IsoChargeSet(int charge);
	virtual ~IsoChargeSet();
	void setSet(vector<DataSet> *set);
	const double * const getNext(int& ,int& );
	int const getLabel(int *setPos);
	int const inline getSize() {return n_points;}
	int const inline getSubSetSize(int ix) {return (*pSet)[ix].getSize();}
	DataSet * inline getSubSet(int ix) {return &((*pSet)[ix]);}
protected:
    int charge;
    vector<DataSet> *pSet;
    int n_points;
    Normalizer norm;
};

#endif /*ISOCHARGESET_H_*/
