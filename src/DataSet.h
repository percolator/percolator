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

#include "Scores.h"
#include "ResultHolder.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include "DescriptionOfCorrect.h"
using namespace std;

// Optional columns in tab delimited input
enum OptionalField {
  SCANNR, EXPMASS, CALCMASS
};

class DataSet {
  public:
    DataSet();
    virtual ~DataSet();
    
    void initFeatureTables(const unsigned int numFeatures);
    
    void inline setLabel(int l) { label_ = l; }
    int inline getLabel() const { return label_; }
    
    void inline setSize(int n) { numSpectra_ = n; }
    int inline getSize() const { return numSpectra_; }
    
    static inline void setCalcDoc(bool on) { calcDOC_ = on; }
    static inline bool getCalcDoc() { return calcDOC_; }
    
    static FeatureNames& getFeatureNames() { return featureNames_; }
    static unsigned getNumFeatures() { return featureNames_.getNumFeatures(); }
    
    void setRetentionTime(map<int, double>& scan2rt) { 
      PSMDescription::setRetentionTime(psms_, scan2rt);
    }
    
    bool writeTabData(ofstream& out);
    
    void print_10features();
    void print_features();
    void print(Scores& test, vector<ResultHolder> & outList);
    
    void fillFeatures(vector<ScoreHolder> &scores);
    void fillFeatures(vector<double*> &features);
    void fillRtFeatures(vector<double*> &rtFeatures);
    
    static double isEnz(const char n, const char c);
    static unsigned int peptideLength(const string& pep);
    static unsigned int cntPTMs(const string& pep);
    // static double isPngasef(const string& peptide, bool isDecoy );

    void readPsm(const std::string line, const unsigned int lineNr,
                 const std::vector<OptionalField>& optionalFields);
    void registerPsm(PSMDescription * myPsm);

  protected:   
    // double isPngasef(const string& peptide);
    static bool calcDOC_;
    const static std::string aaAlphabet_;
    static std::string ptmAlphabet_;
    
    std::vector<PSMDescription*> psms_;
    int label_;
    int numSpectra_;
    static FeatureNames featureNames_;
};

#endif /*DATASET_H_*/
