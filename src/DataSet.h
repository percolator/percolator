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
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <cerrno>

#include "Scores.h"
#include "ResultHolder.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "PSMDescriptionDOC.h"
#include "FeatureNames.h"
#include "DescriptionOfCorrect.h"
#include "FeatureMemoryPool.h"
#include "ProteinProbEstimator.h"

// using char pointers is much faster than istringstream
class TabReader {
 public:
  TabReader(const std::string& line) : f_(line.c_str()), err(0) {
    errno = 0;
  }
  
  void advance(const char* next) {
    if (*next != '\0') {
      f_ = next + 1; // eats up the tab
    } else {
      f_ = next; // prevents pointing over the null byte
    }
  }
  
  void skip() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
    } else {
      advance(pch);
    }
  }
  
  void skip(size_t numSkip) {
    for (size_t i = 0; i < numSkip; ++i) skip();
  }
  
  double readDouble() {
    char* next = NULL;
    errno = 0;
    double d = strtod(f_, &next);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || ((d == HUGE_VAL || d == -HUGE_VAL) && errno == ERANGE)) {
      err = errno ? errno : 1;
    }
    advance(next);
    return d;
  }
  
  int readInt() {
    char* next = NULL;
    errno=0;
    long val = strtol(f_, &next, 10);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || val < INT_MIN || val > INT_MAX) {
      err = errno ? errno : 1;
    }
    advance(next);
    return static_cast<int>(val);
  }
  
  std::string readString() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
      return std::string(f_);
    } else {
      std::string s(f_, static_cast<std::basic_string<char>::size_type>(pch - f_));
      advance(pch);
      return s;
    }
  }

  bool error() { return err != 0; }
 private:
  const char* f_;
  int err;
};

class TabFileValidator {
  public:
    static bool isTabFile(std::string file_name);
    static bool isTabFiles(std::vector<std::string> files);
    std::string getDecoyPrefix(std::vector<std::string> fileList);
    std::string detect_decoy_prefix(std::string file_name);
    bool validateTabFiles(std::vector<std::string> files, std::string* decoy_prefix);
    void getProteinIndex(std::string file_name, int* proteinIndex,int* labelIndex);
    std::string findDecoyPrefix(std::string file_name, int proteinIndex, int labelIndex);
};


// Optional columns in tab delimited input
enum OptionalField {
  SCANNR, EXPMASS, CALCMASS, RETTIME, DELTAMASS
};

class DataSet {
 public:
  DataSet();
  virtual ~DataSet();
  
  void initFeatureTables(const unsigned int numFeatures);
  
  void inline setLabel(int l) { label_ = l; }
  int inline getLabel() const { return label_; }
  
  unsigned int inline getSize() const { return static_cast<unsigned int>(psms_.size()); }
  
  static inline void setCalcDoc(bool on) { calcDOC_ = on; }
  static inline bool getCalcDoc() { return calcDOC_; }
  
  static FeatureNames& getFeatureNames() { return featureNames_; }
  static void resetFeatureNames() { 
    featureNames_ = FeatureNames();
    FeatureNames::resetNumFeatures();
  }
  static unsigned getNumFeatures() { return static_cast<unsigned>(featureNames_.getNumFeatures()); }
  
  bool writeTabData(std::ofstream& out);
  
  void print_10features();
  void print_features();

  void fillScores(std::vector<ScoreHolder>& scores);
  void fillFeatures(std::vector<double*>& features);
  void fillDOCFeatures(std::vector<double*>& features);
  void fillRtFeatures(std::vector<double*>& rtFeatures);
  
  void readPsm(const std::string& line, const unsigned int lineNr,
               const std::vector<OptionalField>& optionalFields, 
               FeatureMemoryPool& featurePool, std::string decoyPrefix);
  static int readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, std::string decoyPrefix);
  
  void registerPsm(PSMDescription* myPsm);
  
 protected:   
  static bool calcDOC_;
  
  std::vector<PSMDescription*> psms_;
  int label_;
  static FeatureNames featureNames_;
};

#endif /*DATASET_H_*/
