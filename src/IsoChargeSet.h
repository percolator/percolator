#ifndef ISOCHARGESET_H_
#define ISOCHARGESET_H_
class IsoChargeSet
{
protected:
    int charge;
    vector<DataSet> *pSet;
    int n_points;
    Normalizer * norm;
public:
	IsoChargeSet(int charge);
	virtual ~IsoChargeSet();
	void setSet(vector<DataSet> *set);
	const double * getNext(int& ,int& );
	int const getLabel(int *setPos);
	inline int const getSize() {return n_points;}
	inline int const getCharge() {return charge;}
	inline int const getSubSetSize(int ix) {return (*pSet)[ix].getSize();}
	inline DataSet * getSubSet(int ix) {return &((*pSet)[ix]);}
	inline Normalizer * getNormalizer() {return norm;} 
};

#endif /*ISOCHARGESET_H_*/
