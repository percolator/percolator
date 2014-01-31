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
#include <boost/foreach.hpp>
#include "ResultHolder.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "SanityCheck.h"
using namespace std;

class SetHandler {
  
  protected:
    vector<DataSet*> subsets;
    int n_examples;
    
    unsigned int getSubsetIndexFromLabel(int label);
    
  public:
    
    SetHandler();
    virtual ~SetHandler();

    void push_back_dataset( DataSet * ds );

    /*void modifyFile(const string& fn, vector<DataSet*> & sets, Scores& sc,
               const string& greet, bool dtaSelect);*/

    void fillTestSet(SetHandler& trainSet, const string& shuffled2FN = "");
    void createXvalSets(vector<SetHandler>& train,
                        vector<SetHandler>& test,
                        const unsigned int xval_fold);
       
    //const double* getFeatures(const int setPos, const int ixPos) const;  
    int readTab(const string& dataFN, SanityCheck *& pCheck);
    void writeTab(const string& dataFN, SanityCheck * pCheck);
    void print(Scores& test, int label, ostream& myout = cout);
    void fillFeatures(vector<ScoreHolder> &scores, int label);
    void fillFeaturesPeptide(vector<ScoreHolder> &scores, int label);

    int const getLabel(int setPos);
    inline int const getSize() { return n_examples; }
    inline int getSizeFromLabel(int label) {
      return (subsets[getSubsetIndexFromLabel(label)]->getSize());
    }
    
    vector<DataSet*> & getSubsets() { return subsets; }
    inline DataSet* getSubset(unsigned int ix) { return (subsets[ix]); }
    inline DataSet* getSubsetFromLabel(int label) {
      return (subsets[getSubsetIndexFromLabel(label)]);
    }
};

#endif /*SETHANDLER_H_*/
