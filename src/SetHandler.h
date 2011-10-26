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
#ifndef SETHANDLER_H_
#define SETHANDLER_H_
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <string>
#include "ResultHolder.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
// #include <rpc/types.h>
// #include <rpc/xdr.h>
// #if defined __WIN32__ or defined __MINGW__
//  #include <xdr_api_mingw.h>
// #endif
#include "percolator_in.hxx"
using namespace std;

namespace percolatorInNs { 
  class target_decoy;
}

class SetHandler {
  protected:
    vector<DataSet*> subsets;
    vector<const double*> examples;
    double* labels;
    double* c_vec;
    int n_examples;
  public:
    SetHandler();
    virtual ~SetHandler();
    void filelessSetup(const unsigned int numFeatures,
                       const unsigned int numSpectra, const int label);



    void push_back_dataset( DataSet * ds );

    void static
    modifyFile(const string& fn, vector<DataSet*> & sets, Scores& sc,
               const string& greet, bool dtaSelect);
    void generateTrainingSet(const double fdr, const double cpos,
                             const double cneg, Scores& sc);
    void setSet();
    void fillTestSet(SetHandler& trainSet, const string& shuffled2FN = "");
    void createXvalSets(vector<SetHandler>& train,
                        vector<SetHandler>& test,
                        const unsigned int xval_fold);
    PSMDescription* getNext(int&, int&);
    const double* getFeatures(const int setPos, const int ixPos) const;
    void readTab(const string& dataFN, const int label);
    void static writeTab(const string& dataFN, const SetHandler& norm,
                         const SetHandler& shuff);
    int const getLabel(int setPos);
    inline int const getTrainingSetSize() {
      return examples.size();
    }
    void print(Scores& test, ostream& myout = cout);
    inline int const getSize() {
      return n_examples;
    }
    inline DataSet* getSubSet(int ix) {
      return (subsets[ix]);
    }
    vector<DataSet*> & getSubsets() {
      return subsets;
    }
    vector<const double*> * getTrainingSet() {
      return &examples;
    }
    inline const double* getLabels() {
      return labels;
    }
    inline const double* getC() {
      return c_vec;
    }
    class Iterator {
      public:
        Iterator(SetHandler* s) :
          sh(s), set(0), ix(-1) {
          ;
        }
        PSMDescription* getNext() {
          return sh->getNext(set, ix);
        }
      private:
        SetHandler* sh;
        int set;
        int ix;
    };
    Iterator getIterator() {
      return Iterator(this);
    }
};

#endif /*SETHANDLER_H_*/
