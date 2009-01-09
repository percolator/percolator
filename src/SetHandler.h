/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-9 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Sciences at the University of Washington. 
 *
 * $Id: SetHandler.h,v 1.28 2009/01/09 14:40:59 lukall Exp $
 *******************************************************************************/
#ifndef SETHANDLER_H_
#define SETHANDLER_H_
#include "PSMDescription.h"

class SetHandler
{
protected:
    vector<DataSet *> subsets;
    vector<const double *> examples;
    double * labels;
    double * c_vec;
    int n_examples;
public:
	SetHandler();
	virtual ~SetHandler();
    void filelessSetup(const unsigned int numFeatures, const unsigned int numSpectra,const int label);
    void readFile(const string & p_fn, const int label);
    void readFile(const string & fn, const string & wc,const bool match);
    void static readFile(const string fn, const int label, vector<DataSet *> & sets,const string & wild = "", const bool match=false);
    void static modifyFile(const string& fn, vector<DataSet *> & sets, Scores& sc , const string& greet, bool dtaSelect);
    void modifyFile(const string& fn, Scores& sc , const string& greet, bool dtaSelect);
    void generateTrainingSet(const double fdr,const double cpos, const double cneg, Scores & sc);
	void setSet();
	void fillTestSet(SetHandler& trainSet,const string& shuffled2FN="");
    void createXvalSets(vector<SetHandler>& train,vector<SetHandler>& test, const unsigned int xval_fold);
	PSMDescription * getNext(int& ,int& );
    const double * getFeatures(const int setPos,const int ixPos) const;
    void readTab(const string & dataFN, const int label);
    void static writeTab(const string &dataFN, const SetHandler& norm, const SetHandler& shuff);
    void readGist(const string & dataFN, const string & labelFN, const int label);
    void static gistWrite(const string & fileNameTrunk,const SetHandler& norm,const SetHandler& shuff);
	int const getLabel(int setPos);
    inline int const getTrainingSetSize() {return examples.size();}
    void print(Scores &test, ostream& myout=cout);    
    inline int const getSize() {return n_examples;}
    inline DataSet * getSubSet(int ix) {return (subsets[ix]);}
    vector<DataSet *> & getSubsets() {return subsets;}
    vector<const double *> * getTrainingSet() {return & examples;}
    inline const double * getLabels() {return labels;}
    inline const double * getC() {return c_vec;}
    class Iterator {
    public:
      Iterator(SetHandler * s) : sh(s), set(0), ix(-1) {;}
      PSMDescription* getNext() {return sh->getNext(set,ix);}
    private:
      SetHandler * sh;
      int set;
      int ix;
    };
    Iterator getIterator() {return Iterator(this);}
};

#endif /*SETHANDLER_H_*/
