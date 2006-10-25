#ifndef CALLER_H_
#define CALLER_H_
class Caller
{
public:
	Caller();
	virtual ~Caller();
	void step(double *w,SetHandler & train);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void readFile(const string fn, const int label, vector<DataSet *> & sets);
    void modifyFile(const string fn, vector<DataSet *> & sets, Scores &sc , const string greet);
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
    double fdr;
    double Cpos;
    double Cneg;
    int niter;
    time_t startTime;
    clock_t startClock;
};

#endif /*CALLER_H_*/
