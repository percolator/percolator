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
#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <cctype>
#include <locale>
#include <queue>
#include <climits>

#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "SanityCheck.h"
#include "PseudoRandom.h"
#include "FeatureMemoryPool.h"

struct PSMDescriptionPriority {
  PSMDescription* psm;
  size_t priority;
  LabelType label;
  
  bool operator<(const PSMDescriptionPriority& psmp) const {
    return (priority < psmp.priority);
  }
};

/*
* ScanId is a tuple of (ScanNr, ExpMass).
*/
typedef std::pair<int, double> ScanId;

struct ScanIdHash {
  // noexcept lets the hash table skip storing the hash code in every node
  size_t operator()(const ScanId& scanId) const noexcept {
    return std::hash<int>()(scanId.first) ^
           (std::hash<double>()(scanId.second) << 1);
  }
};

/*
* SetHandler is a class that provides functionality to handle training,
* testing, Xval data sets, reads/writes from/to a file, prints them.
*
*/
class SetHandler {
 public:
  SetHandler(unsigned int maxPSMs);
  virtual ~SetHandler();
  SetHandler(const SetHandler&) = delete;
  SetHandler& operator=(const SetHandler&) = delete;

  void push_back_dataset(std::unique_ptr<DataSet> ds);

  void setDecoyPrefix(std::string& decoyPrefix) { decoyPrefix_ = decoyPrefix; }

     
  //const double* getFeatures(const int setPos, const int ixPos) const; 
  size_t getMaxPSMs() { return maxPSMs_; }
  // number of PSM lines in the input, counted by the last full read
  size_t getNumPsmLines() const { return numPsmLines_; }
  
  // Reads in tab delimited stream and returns a SanityCheck object based on
  // the presence of default weights. Returns 0 on error, 1 on success.
  int readTab(std::istream& dataStream, SanityCheck& check);
  int readAndScoreTab(std::istream& dataStream, 
    std::vector<double>& rawWeights, Scores& allScores, SanityCheck& check);
  void addQueueToSets(std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet);
  
  void writeTab(const std::string& dataFN, const SanityCheck& check);
  void populateScoresWithPSMs(std::vector<ScoreHolder> &scores, LabelType label);
  void normalizeFeatures(Normalizer& norm);
  
  LabelType const getLabel(int setPos);
  inline int getSizeFromLabel(LabelType label) {
    return static_cast<int>(subsets_[getSubsetIndexFromLabel(label)]->getSize());
  }
  
  inline DataSet* getSubset(unsigned int ix) { return subsets_[ix].get(); }
  inline DataSet* getSubsetFromLabel(LabelType label) {
    return subsets_[getSubsetIndexFromLabel(label)].get();
  }
  
  FeatureMemoryPool& getFeaturePool() { return featurePool_; }
  
  void reset();
  void releaseFeatures();

 protected:
  size_t maxPSMs_;
  size_t numPsmLines_ = 0;
  PsmStore psmStore_;  // owns the PSMs of subsets_ and of -N pass 2
  std::vector<std::unique_ptr<DataSet> > subsets_;
  FeatureMemoryPool featurePool_;
  std::string decoyPrefix_; // Used to determine if a psm is a decoy
  
  unsigned int getSubsetIndexFromLabel(LabelType label);
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
    
  void readPSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, bool& separateSearches,
    std::vector<OptionalField>& optionalFields);
  void readAndScorePSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields, 
    std::vector<double>& rawWeights, Scores& allScores);
};

#endif /*SETHANDLER_H_*/
