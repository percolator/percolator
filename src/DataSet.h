/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: DataSet.h,v 1.46 2008/04/02 00:06:57 lukall Exp $
 *******************************************************************************/
#ifndef DATASET_H_
#define DATASET_H_
#include <string>
#include <set>
#include <vector>
#include <iostream>
using namespace std;
class Scores;
class Normalizer;
class ResultHolder;
class IntraSetRelation;

typedef enum {NO_ENZYME,TRYPSIN,CHYMOTRYPSIN,ELASTASE} Enzyme;


class DataSet
{
 public:
	DataSet();
	virtual ~DataSet();
    static string getFeatureNames();
    static void setFeatureNames(string fn){DataSet::featureNames=fn;}
    void inline setLabel(int l) {label=l;}
    void computeIntraSetFeatures();
    void computeIntraSetFeatures(double *feat,string &pep,set<string> &prots);
    void computeAAFrequencies(const string& pep, double *feat);
    void readSQT(const string fname,IntraSetRelation * intrarel,const string & wild="", bool match=false);
    void modifySQT(const string & outFN, const double *w, Scores * pSc ,const string greet, bool dtaSelect);
    void filelessSetup(const unsigned int numFeatures, const unsigned int numSpectra);
    static inline int getNumFeatures() { return numFeatures; }
    static void setQuadraticFeatures(bool on) { calcQuadraticFeatures=on; }
    static void setCalcIntraSetFeatures(bool on) { calcIntraSetFeatures=on; }
    static void setEnzyme(Enzyme enz) { enzyme=enz; }
    static void setAAFreqencies(bool on) { calcAAFrequencies=on; }
    static void setPTMfeature(bool on) { calcPTMs=on; }      
    static void setNumFeatures();
    static void inline setHitsPerSpectrum(int hits) {hitsPerSpectrum=hits;}
    static inline int rowIx(int row) { return row*numFeatures; }
    double * getFeature() {return feature;}
    const double * getFeatures(const int pos) const;
    int inline getSize() const {return numSpectra;}
    int inline const getLabel() const {return label;}
    double * getNext(int& pos);
    bool writeTabData(ofstream & out, const string & lab);
    void readTabData(ifstream & dataStream, const vector<unsigned int> &ixs);
    bool getGistDataRow(int& pos,string & out);
    void readGistData(ifstream & is, const vector<unsigned int> &ixs);
    void print_10features();
    void print_features();
    void print(Scores& test, vector<ResultHolder> & outList);
protected:
    void readFeatures(const string &in,double *feat,int match,set<string> & proteins, string & pep,bool getIntra);
    string modifyRec(const string record, const set<int>& theMs, const double *w, Scores * pSc, bool dtaSelect);
    static unsigned int peptideLength(const string& pep);
    static unsigned int cntPTMs(const string& pep);
    static unsigned int cntEnz(const string& peptide);
    static double isTryptic(const char n,const char c);
    static double isChymoTryptic(const char n,const char c);
    static double isElastasic(const char n,const char c);
    static double isEnz(const char n,const char c); 
    vector<string> ids;
    static bool calcQuadraticFeatures;
    static bool calcAAFrequencies;
    static Enzyme enzyme;
    static bool calcIntraSetFeatures;
    static bool calcPTMs;
    static int numFeatures;
    static int numRealFeatures;
    static int hitsPerSpectrum;
    static string aaAlphabet;
    static string ptmAlphabet;
    static string featureNames;
    const static int maxNumRealFeatures = 16 + 3 + 20 + 1 + 1; // Normal + Amino acid + PTM + hitsPerSpectrum
    vector<set<string> > proteinIds;
    vector<string> pepSeq;
    int label;
    double *feature;
    int numSpectra;
    string sqtFN;
    string pattern;
    bool doPattern;
    bool matchPattern;
    IntraSetRelation * intra;
};

#endif /*DATASET_H_*/
