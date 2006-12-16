#ifndef CALLER_H_
#define CALLER_H_

typedef enum {NO_XV=0, EACH_STEP, WHOLE} XvType;

class Caller
{
public:
	Caller();
	virtual ~Caller();
    void step(Scores& train,double * w, double Cpos, double Cneg, double fdr);
    void trainEm(double * w);
    void xvalidate_step(double *w);
    void xvalidate(double *w);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, double * weights);
    int run();
protected:
    Normalizer * pNorm;
    AlgIn *svmInput;
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
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    int niter;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    double test_fdr;
    XvType xv_type; 
    vector<Scores> xv_train,xv_test;
    vector<double> xv_fdrs,xv_cposs,xv_cfracs;
    SetHandler normal,shuffled,shuffled2;
    Scores trainset,testset;
};

#endif /*CALLER_H_*/
