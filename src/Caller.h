#ifndef CALLER_H_
#define CALLER_H_


class Caller
{
public:
	Caller();
	virtual ~Caller();
    void step(SetHandler & train,double * w, double Cpos, double Cneg, double fdr);
    void trainEm(SetHandler & set ,double * w, double Cpos, double Cneg, double fdr);
    void xvalidate(vector<DataSet *> &forward,vector<DataSet *> &shuffled, double *w);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void readFile(const string fn, const int label, vector<DataSet *> & sets, IntraSetRelation * intra, bool calc=true);
    void modifyFile(const string fn, vector<DataSet *> & sets, Scores &sc , const string greet);
    void printWeights(ostream & weightStream, double * weights);
    int run();
protected:
    Normalizer * pNorm;
    Scores scores;
    string modifiedFN;
    string forwardFN;
    string shuffledFN;
    string shuffled2FN;
    string rocFN;
    string gistFN;
    string weightFN;
    string call;
    double selectedfdr;
    double selectedCpos;
    double selectedCneg;
    int niter;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold = 3;
    const static double test_fdr = 0.01;
};

#endif /*CALLER_H_*/
