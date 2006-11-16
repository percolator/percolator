#ifndef SETHANDLER_H_
#define SETHANDLER_H_
class SetHandler
{
protected:
    vector<VirtualSet *> subsets;
    vector<const double *> examples;
    double * labels;
    double * c_vec;
    int n_examples;
    int n_pos;
    int n_neg;
    IntraSetRelation * n_intra;
    IntraSetRelation * s_intra;
public:
	SetHandler();
	virtual ~SetHandler();
    void readFile(const string & p_fn,const string & n_fn);
    void readOneFile(const string & fn, const string & wc);
    void static readFile(const string fn, const int label, vector<VirtualSet *> & sets, IntraSetRelation * intra,const string & wild = "", const bool match=false, bool calc=true);
    void static modifyFile(const string& fn, vector<VirtualSet *> & sets, double * w, Scores& sc , const string& greet);
    void modifyFile(const string& forw_fn, const string& shuff_fn, double * w, Scores& sc , const string& greet);
    void generateTrainingSet(const double fdr,const double cpos, const double cneg, const Scores & sc);
	void setSet(vector<VirtualSet *> & pos, vector<VirtualSet *> & neg);
	void setSet(vector<VirtualSet *> & sets);
	void setSet();
	void fillTestSet(SetHandler& trainSet,const string& shuffled2FN="");
    void createXvalSets(vector<SetHandler>& train,vector<SetHandler>& test, const unsigned int xval_fold);
	const double * getNext(int& ,int& );
    const double * getFeatures(const int setPos,const int ixPos);
    void readGist(const string & dataFN, const string & labelFN);
    void gistWrite(const string & fileNameTrunk);
	int const getLabel(int setPos);
    inline int const getTrainingSetSize() {return examples.size();}
    
    inline int const getSize() {return n_examples;}
//    inline int const getPositiveSize() {return n_pos;}
    inline int const getNegativeSize() {return n_neg;}
    inline VirtualSet * getSubSet(int ix) {return (subsets[ix]);}
    vector<VirtualSet *> & getSubsets() {return subsets;}
    vector<const double *> * getTrainingSet() {return & examples;}
    inline const double * getLabels() {return labels;}
    inline const double * getC() {return c_vec;}

};

#endif /*SETHANDLER_H_*/
