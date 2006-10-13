#ifndef SETHANDLER_H_
#define SETHANDLER_H_
class SetHandler
{
protected:
//    int charge;
    vector<DataSet *> subsets;
    int n_points;
    Normalizer * norm;
public:
	SetHandler();
	virtual ~SetHandler();
	void setSet(DataSet &pos, DataSet &neg);
	const double * getNext(int& ,int& );
    void gistWrite(const string & fileNameTrunk);
	int const getLabel(int *setPos);
	inline int const getSize() {return n_points;}
//	inline int const getCharge() {return charge;}
	inline int const getSubSetSize(int ix) {return subsets[ix]->getSize();}
	inline DataSet * getSubSet(int ix) {return (subsets[ix]);}
	inline Normalizer * getNormalizer() {return norm;} 
};

#endif /*SETHANDLER_H_*/
