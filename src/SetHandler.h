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

#include "ResultHolder.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "SanityCheck.h"

using namespace std;

/*
* SetHandler is a class that provides functionality to handle training,
* testing, Xval data sets, reads/writes from/to a file, prints them.
*
*/
class SetHandler {    
 public:
  SetHandler();
  virtual ~SetHandler();

  void push_back_dataset(DataSet* ds);
     
  //const double* getFeatures(const int setPos, const int ixPos) const; 
  
  // Reads in tab delimited stream and returns a SanityCheck object based on
  // the presence of default weights. Returns 0 on error, 1 on success.
  int readTab(istream& dataStream, SanityCheck*& pCheck);
  void writeTab(const string& dataFN, SanityCheck* pCheck);
  void print(Scores& test, int label, ostream& myout = cout);
  void fillFeatures(vector<ScoreHolder> &scores, int label);

  int const getLabel(int setPos);
  inline int getSizeFromLabel(int label) {
    return (subsets_[getSubsetIndexFromLabel(label)]->getSize());
  }
  
  vector<DataSet*>& getSubsets() { return subsets_; }
  inline DataSet* getSubset(unsigned int ix) { return (subsets_[ix]); }
  inline DataSet* getSubsetFromLabel(int label) {
    return (subsets_[getSubsetIndexFromLabel(label)]);
  }
  
 protected:
  vector<DataSet*> subsets_;
  
  unsigned int getSubsetIndexFromLabel(int label);
  static inline std::string &rtrim(std::string &s);
};

#endif /*SETHANDLER_H_*/
