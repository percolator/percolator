#ifndef CALLER_H_
#define CALLER_H_
class Caller
{
public:
	Caller();
	virtual ~Caller();
	inline void setSet(SetHandler *p){pSet=p;}
	void step(double *);
    bool parseOptions(int argc, char **argv);
    int run();
protected:
    SetHandler *pSet;
    Scores scores;
    string forwardFN;
    string shuffledFN;
    string shuffled2FN;
    string rocFN;
    string gistFN;
    string weightFN;
    double fdr;
    int nitter;
};

#endif /*CALLER_H_*/
