/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <boost/foreach.hpp>
#include "Scores.h"
#include "ResultHolder.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include "DescriptionOfCorrect.h"
using namespace std;

class DataSet {
  public:
    DataSet();
    virtual ~DataSet();
    
    void initFeatureTables(const unsigned int numFeatures);
    
    void inline setLabel(int l) { label = l; }
    int inline getLabel() const { return label; }
    
    void inline setSize(int n) { numSpectra = n; }
    int inline getSize() const { return numSpectra; }
    
    static inline void setCalcDoc(bool on) { calcDOC = on; }
    static inline bool getCalcDoc() { return calcDOC; }
    
    static void setIsotopeMass(bool on) { isotopeMass = on; }
    
    static FeatureNames& getFeatureNames() { return featureNames; }
    static unsigned getNumFeatures() { return featureNames.getNumFeatures(); }
    
    void setRetentionTime(map<int, double>& scan2rt) { PSMDescription::setRetentionTime(psms, scan2rt); }
    
    bool writeTabData(ofstream& out);
    
    void print_10features();
    void print_features();
    void print(Scores& test, vector<ResultHolder> & outList);
    
    void fillFeatures(vector<ScoreHolder> &scores);
    void fillFeaturesPeptide(vector<ScoreHolder> &scores);
    void fillFeatures(vector<double*> &features);
    void fillRtFeatures(vector<double*> &rtFeatures);
    
    static double isEnz(const char n, const char c);
    static unsigned int peptideLength(const string& pep);
    static unsigned int cntPTMs(const string& pep);
    // static double isPngasef(const string& peptide, bool isDecoy );

    void readPsm(ifstream & dataStream, const std::string line);
    void registerPsm(PSMDescription * myPsm);

  protected:   
    // double isPngasef(const string& peptide);
    static bool calcDOC;
    static bool isotopeMass;
    const static string aaAlphabet;
    static string ptmAlphabet;
    const static int maxNumRealFeatures = 16 + 3 + 20 * 3 + 1 + 1 + 3; // Normal + Amino acid + PTM + hitsPerSpectrum + doc
    
    vector<PSMDescription*> psms;
    int label;
    int numSpectra;
    string fileId;
    static FeatureNames featureNames;
};

#endif /*DATASET_H_*/
