/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

 *******************************************************************************/
#ifndef DATASET_H_
#define DATASET_H_
#include <string>
#include <set>
#include <map>
#include <vector>
#include <iostream>
#include "PSMDescription.h"
#include "FeatureNames.h"
using namespace std;
class Scores;
class Normalizer;
class ResultHolder;

typedef enum {NO_ENZYME,TRYPSIN,CHYMOTRYPSIN,ELASTASE} Enzyme;


class DataSet
{
 public:
	DataSet();
	virtual ~DataSet();
    void inline setLabel(int l) {label=l;}
    void computeAAFrequencies(const string& pep, double *feat);
    void readSQT(const string fname,const string & wild="", bool match=false);
    void modifySQT(const string & outFN, Scores * pSc ,const string greet, bool dtaSelect);
    void initFeatureTables(const unsigned int numFeatures, const unsigned int numSpectra, bool regresionTable = false);
    static FeatureNames& getFeatureNames() { return featureNames; }
    static void setQuadraticFeatures(bool on) { calcQuadraticFeatures=on; }
    static void setCalcDoc(bool on) { calcDOC=on; }
    static bool getCalcDoc() { return calcDOC; }
    static void setEnzyme(Enzyme enz) { enzyme=enz; }
    static void setAAFreqencies(bool on) { calcAAFrequencies=on; }
    static void setPTMfeature(bool on) { calcPTMs=on; }
    static void setPNGaseF(bool on) { pngasef=on; }
    static void setIsotopeMass(bool on) { isotopeMass=on; }
    static void setNumFeatures(bool doc);
    static void inline setHitsPerSpectrum(int hits) {hitsPerSpectrum=hits;}
    static inline int rowIx(int row) { return row*FeatureNames::getNumFeatures(); }
    double * getFeature() {return feature;}
    const double * getFeatures(const int pos) const;
    int inline getSize() const {return numSpectra;}
    int inline const getLabel() const {return label;}
    PSMDescription* getNext(int& pos);
    void setRetentionTime(map<int,double>& scan2rt) {PSMDescription::setRetentionTime(psms,scan2rt);}
    bool writeTabData(ofstream & out, const string & lab);
    void readTabData(ifstream & dataStream, const vector<unsigned int> &ixs);
    bool getGistDataRow(int& pos,string & out);
    void readGistData(ifstream & is, const vector<unsigned int> &ixs);
    void print_10features();
    void print_features();
    void print(Scores& test, vector<ResultHolder> & outList);
    static double isEnz(const char n,const char c);
protected:
    void readFeatures(const string &in,PSMDescription &psm,int match);
    string modifyRec(const string record, int& row, const set<int>& theMs, Scores * pSc, bool dtaSelect);
    static unsigned int peptideLength(const string& pep);
    static unsigned int cntPTMs(const string& pep);
    static unsigned int cntEnz(const string& peptide);
    static double isTryptic(const char n,const char c);
    static double isChymoTryptic(const char n,const char c);
    static double isElastasic(const char n,const char c);
    double isPngasef(const string& peptide);
    static bool calcQuadraticFeatures;
    static bool calcAAFrequencies;
    static Enzyme enzyme;
    static bool calcPTMs;
    static bool calcDOC;
    static bool isotopeMass;
    static int hitsPerSpectrum;
    static bool pngasef;
    static string reversedFeaturePattern;
    static string aaAlphabet;
    static string ptmAlphabet;
    const static int maxNumRealFeatures = 16 + 3 + 20*3 + 1 + 1 + 3; // Normal + Amino acid + PTM + hitsPerSpectrum + doc
    vector<PSMDescription> psms;
    int label;
    double *feature,*regressionFeature;
    int numSpectra;
    string sqtFN;
    string pattern;
    string fileId;
    bool doPattern;
    bool matchPattern;
    static FeatureNames featureNames;
};

#endif /*DATASET_H_*/
