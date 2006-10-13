#ifndef DATASET_H_
#define DATASET_H_

class DataSet
{
protected:
    int line2fields(char * s, vector<string> *fields);
    double *feature;
    int n_examples;
    vector<string> ids;
    vector<int> charge;
	string delim;
	int label;
    static int numFeatures;
    static bool calcQuadraticFeatures;
    const static int numRealFeatures = 14;
public:
	DataSet();
	virtual ~DataSet();
	double * getFeature() {return feature;}
	int inline getSize() {return n_examples;}
    int inline getLabel() {return label;}
    bool getGistDataRow(int& pos,string & out);
	void inline setLabel(int l) {label=l;}
    double * getNext(int& pos);
    void read_sqt(string & fname);
    void print_features();
    static double isTryptic(const string & str);
    static inline int getNumFeatures() { return numFeatures; }
    static void setQuadraticFeatures(bool on)
      { numFeatures=(on?numRealFeatures*(numRealFeatures+1)/2:numRealFeatures); calcQuadraticFeatures=on;}
    static inline int rowIx(int row) { return row*numFeatures; }
};

#endif /*DATASET_H_*/
