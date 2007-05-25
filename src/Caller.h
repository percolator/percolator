/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Caller.h,v 1.34 2007/05/25 00:09:52 lukall Exp $
 *******************************************************************************/
#ifndef CALLER_H_
#define CALLER_H_

typedef enum {NO_XV=0, EACH_STEP, WHOLE} XvType;

class Caller
{
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
    void readFiles(bool &doSingleFile, bool &separateShuffledTestSetHandler, bool &separateShuffledThresholdSetHandler);
    void filelessSetup(const unsigned int sets, const unsigned int numFeatures, const unsigned int numSpectra);
    void fillFeatureSets(bool &separateShuffledTestSetHandler, bool &separateShuffledThresholdSetHandler);    
    void preIterationSetup();    
    int run();
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
