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

#ifndef SETHANDLER_H_
#define SETHANDLER_H_

#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <cctype>
#include <locale>
#include <queue>
#include <climits>

#include "ResultHolder.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "SanityCheck.h"
#include "PseudoRandom.h"
#include "DescriptionOfCorrect.h"
#include "FeatureMemoryPool.h"

using namespace std;

struct PSMDescriptionPriority {
  PSMDescription* psm;
  size_t priority;
  int label;
  
  bool operator<(const PSMDescriptionPriority& psmp) const {
    return (priority < psmp.priority);
  }
};

typedef std::pair<int, double> ScanId;

/*
* SetHandler is a class that provides functionality to handle training,
* testing, Xval data sets, reads/writes from/to a file, prints them.
*
*/
class SetHandler {    
 public:
  SetHandler(unsigned int maxPSMs);
  virtual ~SetHandler();

  void push_back_dataset(DataSet* ds);
     
  //const double* getFeatures(const int setPos, const int ixPos) const; 
  size_t getMaxPSMs() { return maxPSMs_; }
  
  // Reads in tab delimited stream and returns a SanityCheck object based on
  // the presence of default weights. Returns 0 on error, 1 on success.
  int readTab(istream& dataStream, SanityCheck*& pCheck, std::string decoyPrefix);
  int readAndScoreTab(istream& dataStream, 
    std::vector<double>& rawWeights, Scores& allScores, SanityCheck*& pCheck, std::string decoyPrefix);
  void addQueueToSets(std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet);
  
  void writeTab(const string& dataFN, SanityCheck* pCheck);
  void populateScoresWithPSMs(vector<ScoreHolder> &scores, int label);
  void normalizeFeatures(Normalizer*& pNorm);
  void normalizeDOCFeatures(Normalizer* pNorm);
  
  int const getLabel(int setPos);
  inline int getSizeFromLabel(int label) {
    return static_cast<int>(subsets_[getSubsetIndexFromLabel(label)]->getSize());
  }
  
  inline DataSet* getSubset(unsigned int ix) { return (subsets_[ix]); }
  inline DataSet* getSubsetFromLabel(int label) {
    return (subsets_[getSubsetIndexFromLabel(label)]);
  }
  
  static void deletePSMPointer(PSMDescription* psm);
  
  FeatureMemoryPool& getFeaturePool() { return featurePool_; }
  
  void reset();

  std::string getDecoyPrefix(std::vector<std::string> fileList);
  std::string detect_decoy_prefix(std::string file_name);
  
 protected:
  size_t maxPSMs_;
  vector<DataSet*> subsets_;
  FeatureMemoryPool featurePool_;
  
  unsigned int getSubsetIndexFromLabel(int label);
  static inline std::string &rtrim(std::string &s);
  
  int getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields);
  bool isDefaultDirectionLine(const std::string& defaultDirectionLine);
  int getNumFeatures(const std::string& line, int optionalFieldCount);
  void getFeatureNames(const std::string& headerLine, int numFeatures, 
    int optionalFieldCount, FeatureNames& featureNames);
  bool getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values);
  ScanId getScanId(const std::string& psmLine, int& label,
    std::vector<OptionalField>& optionalFields, unsigned int lineNr);
    
  void readPSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, bool& separateSearches,
    std::vector<OptionalField>& optionalFields, std::string decoyPrefix);
  void readAndScorePSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields, 
    std::vector<double>& rawWeights, Scores& allScores, std::string decoyPrefix);
};

#endif /*SETHANDLER_H_*/
