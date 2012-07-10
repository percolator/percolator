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

#ifdef WIN32
#ifndef isfinite
#define isfinite _finite
#endif
#endif
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <string>
#include "Scores.h"
#include "ResultHolder.h"
#include "MassHandler.h"
#include "Enzyme.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include <boost/foreach.hpp>
#include "percolator_in.hxx"
using namespace std;

using namespace std;
class Scores;
class Normalizer;
class ResultHolder;

namespace percolatorInNs { 
  class target_decoy;
}

class DataSet {
  public:
    DataSet();
    virtual ~DataSet();
    void inline setLabel(int l) {
      label = l;
    }
    void readTargetDecoy(const ::percolatorInNs::target_decoy & td, unsigned int numFeatures );
    
    void initFeatureTables(const unsigned int numFeatures,bool regresionTable = false);
    
    static FeatureNames& getFeatureNames() {
      return featureNames;
    }
    
    static void setCalcDoc(bool on) {
      calcDOC = on;
    }
    static bool getCalcDoc() {
      return calcDOC;
    }
    
    static void setIsotopeMass(bool on) {
      isotopeMass = on;
    }
    
    static void setNumFeatures(bool doc);

    //const double* getFeatures(const int pos) const;
    
    int inline getSize() const {
      return numSpectra;
    }
    
    const int inline  getLabel() const {
      return label;
    }
    
    PSMDescription* getNext(int& pos);
    
    void setRetentionTime(map<int, double>& scan2rt) {
      PSMDescription::setRetentionTime(psms, scan2rt);
    }
    
    bool writeTabData(ofstream& out, const string& lab);
    void readTabData(ifstream& dataStream, const vector<unsigned int> &ixs);
    void print_10features();
    void print_features();
    void print(Scores& test, vector<ResultHolder> & outList);
    static double isEnz(const char n, const char c);
    void readFragSpectrumScans( const ::percolatorInNs::fragSpectrumScan & fss);
    static unsigned int peptideLength(const string& pep);
    static unsigned int cntPTMs(const string& pep);
    static double isPngasef(const string& peptide, bool isDecoy );
    void readPsm(const ::percolatorInNs::peptideSpectrumMatch &psm, unsigned scanNumber );

  protected:
    
    double isPngasef(const string& peptide);
    static bool calcDOC;
    static bool isotopeMass;
    static string reversedFeaturePattern;
    const static string aaAlphabet;
    static string ptmAlphabet;
    const static int maxNumRealFeatures = 16 + 3 + 20 * 3 + 1 + 1 + 3; // Normal + Amino acid + PTM + hitsPerSpectrum + doc
    vector<PSMDescription*> psms;
    int label;
    int numSpectra;
    string sqtFN;
    string pattern;
    string fileId;
    bool doPattern;
    bool matchPattern;
    static FeatureNames featureNames;
    bool regresionTable;
};

#endif /*DATASET_H_*/
