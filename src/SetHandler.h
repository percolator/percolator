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
    IntraSetRelation * intra;
public:
	SetHandler();
	virtual ~SetHandler();
    void readFile(const string & p_fn, const int label);
    void readFile(const string & fn, const string & wc,const bool match);
    void static readFile(const string fn, const int label, vector<DataSet *> & sets, IntraSetRelation * intra,const string & wild = "", const bool match=false, bool calc=true);
    void static modifyFile(const string& fn, vector<DataSet *> & sets, double * w, Scores& sc , const string& greet, bool dtaSelect);
    void modifyFile(const string& fn, double * w, Scores& sc , const string& greet, bool dtaSelect);
    void generateTrainingSet(const double fdr,const double cpos, const double cneg, const Scores & sc);
	void setSet();
	void fillTestSet(SetHandler& trainSet,const string& shuffled2FN="");
    void createXvalSets(vector<SetHandler>& train,vector<SetHandler>& test, const unsigned int xval_fold);
	const double * getNext(int& ,int& ) const;
    const double * getFeatures(const int setPos,const int ixPos) const;
    void readGist(const string & dataFN, const string & labelFN, const int label);
    void static gistWrite(const string & fileNameTrunk,const SetHandler& norm,const SetHandler& shuff,const SetHandler& shuff2);
	int const getLabel(int setPos);
    inline int const getTrainingSetSize() {return examples.size();}
    void print(Scores &test);    
    inline int const getSize() {return n_examples;}
    inline DataSet * getSubSet(int ix) {return (subsets[ix]);}
    vector<DataSet *> & getSubsets() {return subsets;}
    vector<const double *> * getTrainingSet() {return & examples;}
    inline const double * getLabels() {return labels;}
    inline const double * getC() {return c_vec;}

};

#endif /*SETHANDLER_H_*/
