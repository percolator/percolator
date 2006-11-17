#ifndef VIRTUALSET_H_
#define VIRTUALSET_H_

class IntraSetRelation;

class VirtualSet
{
public:
	VirtualSet(const VirtualSet& mother,const int fold,const int ix);
	VirtualSet();
	virtual ~VirtualSet();

	double * getFeature() {return feature;}
    const double * getFeatures(const int pos);
	int inline getSize() {return n_examples;}
    int inline const getLabel() {return label;}
    double * getNext(int& pos);
    bool getGistDataRow(int& pos,string & out);
    void print_10features();
    void print_features();
    inline bool isNormalized() {return *normalizedFlag;}
    inline void setNormalized() {*normalizedFlag=true;}
    
protected:
    int label;
    double *feature;
    int n_examples;
    string sqtFN;
    IntraSetRelation * intra;
    bool *normalizedFlag;
};

#endif /*VIRTUALSET_H_*/
