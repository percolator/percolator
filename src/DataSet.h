#ifndef DATASET_H_
#define DATASET_H_

class Scores;
class Normalizer;
class IntraSetRelation;

class DataSet
{
protected:
    int label;
    double *feature;
    int n_examples;
    int line2fields(string & s, vector<string>  * fields);
    vector<string> ids;
    static bool calcQuadraticFeatures;
    static bool calcTrypticFeatures;
    static bool chymoInsteadOfTryptic;
    static bool calcIntraSetFeatures;
    static int numFeatures;
    static int numRealFeatures;
    static string featureNames;
    const static int maxNumRealFeatures = 16;
    string sqtFN;
    IntraSetRelation * intra;
    vector<set<string> > proteinIds;
    vector<string> pepSeq;
public:
	DataSet();
	virtual ~DataSet();
	double * getFeature() {return feature;}
    const double * getFeatures(const int pos);
	int inline getSize() {return n_examples;}
    int inline getLabel() {return label;}
    static string getFeatureNames();
    static void setFeatureNames(string fn){DataSet::featureNames=fn;}
    bool getGistDataRow(int& pos,string & out);
    void readGistData(ifstream & is, vector<unsigned int> ixs);
	void inline setLabel(int l) {label=l;}
    double * getNext(int& pos);
    void computeIntraSetFeatures();
    void computeIntraSetFeatures(double *feat,string &pep,set<string> &prots);
    void readFeatures(const string &in,double *feat,int match,set<string> & proteins, string & pep,bool getIntra);
    void read_sqt(const string fname,IntraSetRelation * intrarel, const string & wild="", bool match=false);
    string modifyRec(const string record, int mLines, const double *w, Scores * pSc);
    void modify_sqt(const string & outFN, const double *w, Scores * pSc ,const string greet);
    void print_10features();
    void print_features();
    static double isTryptic(const string & str);
    static double isChymoTryptic(const string & str);
    static double isEnz(const string & str) {return (chymoInsteadOfTryptic?
                                                     isChymoTryptic(str):
                                                     isTryptic(str));}
    static inline int getNumFeatures() { return numFeatures; }
    static void setQuadraticFeatures(bool on)
      { calcQuadraticFeatures=on;}
    static void setTrypticFeatures(bool on)
      { calcTrypticFeatures=on;}
    static void setChymoTrypticFeatures(bool on)
      { chymoInsteadOfTryptic=on;}
    static void setNumFeatures();
    static inline int rowIx(int row) { return row*numFeatures; }
};

#endif /*DATASET_H_*/
