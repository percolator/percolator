#ifndef DATASET_H_
#define DATASET_H_
#include <vector>
#include <string>
using namespace std;

#define NUM_FEATURE 8
#define ROW(x) x*NUM_FEATURE

class DataSet
{
protected:
//    int line(string s) {
    int line2fields(char * s, vector<string> *fields);
    double *feature;
    int n_feature;
    vector<string> ids;
    vector<int> charge;
	string delim;
	int label;
public:
	DataSet();
	virtual ~DataSet();
	double * getFeature() {return feature;}
	int inline getSize() {return n_feature;}
	int inline getLabel() {return label;}
	void inline setLabel(int l) {label=l;}
//    vector<int> * getCharges() {return &charge;}
    int getIsoChargeSize(int c);
    double * getNext(const int charge,int& pos);
    void read_sqt(char* fname);
    void print_features();
};

#endif /*DATASET_H_*/
