#ifndef SETHANDLER_H_
#define SETHANDLER_H_
class SetHandler
{
protected:
    vector<DataSet *> subsets;
    vector<const double *> examples;
    double * labels;
    double * c_vec;
    int n_examples;
    int n_pos;
    int n_neg;
public:
	SetHandler();
	virtual ~SetHandler();
    void static readFile(const string fn, const int label, vector<DataSet *> & sets, IntraSetRelation * intra, bool calc=true);
    void static modifyFile(const string fn, vector<DataSet *> & sets, double * w, Scores &sc , const string greet);
    void generateTrainingSet(const double fdr,const double cpos, const double cneg, const Scores & sc);
	void setSet(vector<DataSet *> & pos, vector<DataSet *> & neg);
	const double * getNext(int& ,int& );
    const double * getFeatures(const int setPos,const int ixPos);
    void static readGist(const string dataFN, const string labelFN, vector<DataSet *> & poss, vector<DataSet *> & negs);    void gistWrite(const string & fileNameTrunk);
	int const getLabel(int setPos);
    inline int const getTrainingSetSize() {return examples.size();}
    
    inline int const getSize() {return n_examples;}
//    inline int const getPositiveSize() {return n_pos;}
    inline int const getNegativeSize() {return n_neg;}
    inline DataSet * getSubSet(int ix) {return (subsets[ix]);}
    vector<const double *> * getTrainingSet() {return & examples;}
    inline const double * getLabels() {return labels;}
    inline const double * getC() {return c_vec;}
};

#endif /*SETHANDLER_H_*/
