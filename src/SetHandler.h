#ifndef SETHANDLER_H_
#define SETHANDLER_H_
class SetHandler
{
protected:
//    int charge;
    vector<DataSet *> subsets;
    int n_examples;
    int n_pos;
    int n_neg;
    Normalizer * norm;
public:
	SetHandler();
	virtual ~SetHandler();
	void setSet(vector<DataSet *> & pos, vector<DataSet *> & neg);
	const double * getNext(int& ,int& );
    const double * getFeatures(const int setPos,const int ixPos);
    void gistWrite(const string & fileNameTrunk);
	int const getLabel(int setPos);
    inline int const getSize() {return n_examples;}
    inline int const getPositiveSize() {return n_pos;}
    inline int const getNegativeSize() {return n_neg;}
//	inline int const getCharge() {return charge;}
//	inline int const getSubSetSize(int ix) {return subsets[ix]->getSize();}
	inline DataSet * getSubSet(int ix) {return (subsets[ix]);}
	inline Normalizer * getNormalizer() {return norm;} 
};

#endif /*SETHANDLER_H_*/
