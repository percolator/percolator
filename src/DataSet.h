#ifndef DATASET_H_
#define DATASET_H_
#include <vector>
#include <string>
using namespace std;

class DataSet
{
protected:
//    int line(string s) {
    int line2fields(char * s, vector<string> *fields);
    double **feature;
    int n_feature;
public:
	DataSet();
	virtual ~DataSet();
    void read_sqt(char* fname);
    void print_features();
    const static int ONLY_CHARGE = 2;
	string delim;
};

#endif /*DATASET_H_*/
