/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Caller.h,v 1.49 2008/08/25 14:53:52 lukall Exp $
 *******************************************************************************/
#ifndef CALLER_H_
#define CALLER_H_

#include <time.h>
#include "SanityCheck.h"

class Caller
{
public:
    enum XvType {NO_XV=0, EACH_STEP, WHOLE};
    enum SetHandlerType {NORMAL=0,SHUFFLED,SHUFFLED_TEST,SHUFFLED_THRESHOLD};
public:
	Caller();
	virtual ~Caller();
    void readRetentionTime(string filename);
    void step(Scores& train, vector<double>& w, double Cpos, double Cneg, double fdr);
    void train(vector<vector<double> >& w);
    void trainEm(vector<vector<double> >& w);
    int xv_step(vector<vector<double> >& w);
    void xvalidate(vector<vector<double> >& w);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, vector<double>& w);
    void readWeights(istream & weightStream, vector<double>& w);
    void readFiles(bool &doSingleFile);
    void filelessSetup(const unsigned int numFeatures, const unsigned int numSpectra, char ** fetureNames, double pi0);
    void fillFeatureSets();    
    int preIterationSetup(vector<vector<double> >& w);
    Scores* getFullSet() {return &fullset;}    
    int run();
    SetHandler * getSetHandler(SetHandlerType sh) {
        switch(sh) {
           case NORMAL: return &normal;
           case SHUFFLED: return &shuffled;
           case SHUFFLED_TEST: return NULL;
           case SHUFFLED_THRESHOLD: return NULL;
           default: return NULL;
        }
    }
protected:
    Normalizer * pNorm;
    SanityCheck * pCheck;
    AlgIn *svmInput;
    string modifiedFN;
    string modifiedDecoyFN;
    string forwardFN;
    string decoyFN;
    string decoyWC;
    string rocFN;
    string gistFN;
    string tabFN;
    string weightFN;
    string call;
    string spectrumFile;
    string decoyOut;
    bool gistInput;
    bool tabInput;
    bool dtaSelect;
    bool docFeatures;
    bool reportPerformanceEachIteration;
    double test_fdr;
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    double threshTestRatio;    
    double trainRatio;    
    unsigned int niter;
    unsigned int seed;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    XvType xv_type; 
    vector<Scores> xv_train,xv_test;
    vector<double> xv_cposs,xv_cfracs;
    SetHandler normal,shuffled; //,shuffledTest,shuffledThreshold;
    Scores fullset; //,thresholdset;
    map<int,double> scan2rt; 
};

#endif /*CALLER_H_*/
