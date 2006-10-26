#ifndef DATASET_H_
#define DATASET_H_

class Scores;
class Normalizer;

class DataSet
{
protected:
    int label;
    double *feature;
    int n_examples;
    int line2fields(string & s, vector<string>  * fields);
    vector<string> ids;
    vector<int> charge;
    static bool calcQuadraticFeatures;
    static bool calcTrypticFeatures;
    static bool calcIntraSetFeatures;
    static int numFeatures;
    static int numRealFeatures;
    const static int maxNumRealFeatures = 16;
    string sqtFN;
public:
	DataSet();
	virtual ~DataSet();
	double * getFeature() {return feature;}
    const double * getFeatures(const int pos);
	int inline getSize() {return n_examples;}
    int inline getLabel() {return label;}
    static void getFeatureNames(string & out);
    bool getGistDataRow(int& pos,string & out);
	void inline setLabel(int l) {label=l;}
    double * getNext(int& pos);
    void read_sqt(const string fname);
    void modify_sqt(const string out, vector<double> & sc, vector<double> & fdr, const string greet);
    void print_features();
    static double isTryptic(const string & str);
    static double isChymoTryptic(const string & str);
    static inline int getNumFeatures() { return numFeatures; }
    static void setQuadraticFeatures(bool on)
      { calcQuadraticFeatures=on;}
    static void setTrypticFeatures(bool on)
      { calcTrypticFeatures=on;}
    static void setNumFeatures();
    static inline int rowIx(int row) { return row*numFeatures; }
};

#endif /*DATASET_H_*/
