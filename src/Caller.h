#ifndef CALLER_H_
#define CALLER_H_
class Caller
{
public:
	Caller();
	virtual ~Caller();
	inline void setSet(IsoChargeSet *p){pSet=p;}
	void step(double *);
    bool parseOptions(int argc, char **argv);
    int run();
protected:
    IsoChargeSet *pSet;
    Scores scores;
    string forwardFN;
    string shuffledFN;
    string shuffled2FN;
    string rocFN;
    string gistFN;
    string weightFN;
    bool doRoc;
    double fdr;
};

#endif /*CALLER_H_*/
