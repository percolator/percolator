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
    void printWeights(ostream & weightStream, double * weights);
    int run();
protected:
    Normalizer * pNorm;
    Scores scores;
    string modifiedFN;
    string modifiedShuffledFN;
    string forwardFN;
    string shuffledFN;
    string shuffled2FN;
    string shuffledWC;
    string rocFN;
    string gistFN;
    string weightFN;
    string call;
    bool gistInput;
    double selectedfdr;
    double selectedCpos;
    double selectedCneg;
    int niter;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    const static double test_fdr;
};

#endif /*CALLER_H_*/
