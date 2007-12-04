/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Caller.h,v 1.37 2007/12/04 01:47:52 lukall Exp $
 *******************************************************************************/
#ifndef CALLER_H_
#define CALLER_H_

class Caller
{
public:
    enum XvType {NO_XV=0, EACH_STEP, WHOLE};
    enum SetHandlerType {NORMAL=0,SHUFFLED,SHUFFLED_TEST,SHUFFLED_THRESHOLD};
public:
	Caller();
	virtual ~Caller();
    void step(Scores& train,Scores& thresh,double * w, double Cpos, double Cneg, double fdr);
    void train(double * w);
    void trainEm(double * w);
    void xvalidate_step(double *w);
    void xvalidate(double *w);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, double * weights);
    void readWeights(istream & weightStream, double * weights);
    void readFiles(bool &doSingleFile, bool &separateShuffledTestSetHandler, bool &separateShuffledThresholdSetHandler);
    void filelessSetup(const unsigned int sets, const unsigned int numFeatures, const unsigned int numSpectra, char ** fetureNames, double pi0);
    void fillFeatureSets(bool &separateShuffledTestSetHandler, bool &separateShuffledThresholdSetHandler);    
    int preIterationSetup(double * w);
    Scores* getTestSet() {return &testset;}    
    int run();
    SetHandler * getSetHandler(SetHandlerType sh) {
        switch(sh) {
           case NORMAL: return &normal;
           case SHUFFLED: return &shuffled;
           case SHUFFLED_TEST: return &shuffledTest;
           case SHUFFLED_THRESHOLD: return &shuffledThreshold;
           default: return NULL;
        }
    }
protected:
    Normalizer * pNorm;
    AlgIn *svmInput;
    string modifiedFN;
    string modifiedShuffledFN;
    string forwardFN;
    string shuffledTrainFN;
    string shuffledThresholdFN;
    string shuffledTestFN;
    string shuffledWC;
    string rocFN;
    string gistFN;
    string weightFN;
    string initWeightFN;
    string call;
    bool gistInput;
    bool dtaSelect;
    bool thresholdCalulationOnTrainSet;
    bool reportPerformanceEachIteration;
    double test_fdr;
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    double threshTestRatio;    
    double trainRatio;    
    int niter;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    XvType xv_type; 
    vector<Scores> xv_train,xv_test;
    vector<double> xv_fdrs,xv_cposs,xv_cfracs;
    SetHandler normal,shuffled,shuffledTest,shuffledThreshold;
    Scores trainset,testset,thresholdset;
};

#endif /*CALLER_H_*/
