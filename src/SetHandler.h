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

#ifdef XML_SUPPORT
  #include "percolator_in.hxx"
#endif //XML_SUPPORT
using namespace std;

#ifdef XML_SUPPORT
namespace percolatorInNs { 
  class target_decoy;
}
#endif //XML_SUPPORT

class SetHandler {
  
  protected:
    vector<DataSet*> subsets;
    int n_examples;
    
  public:
    
    SetHandler();
    virtual ~SetHandler();
    void filelessSetup(const unsigned int numFeatures,
                       const unsigned int numSpectra, const int label);

    void push_back_dataset( DataSet * ds );

    /*void modifyFile(const string& fn, vector<DataSet*> & sets, Scores& sc,
               const string& greet, bool dtaSelect);*/
    
    void setSet();
    void fillTestSet(SetHandler& trainSet, const string& shuffled2FN = "");
    void createXvalSets(vector<SetHandler>& train,
                        vector<SetHandler>& test,
                        const unsigned int xval_fold);
    
    PSMDescription* getNext(int&, int&);
    
    //const double* getFeatures(const int setPos, const int ixPos) const;
    
    void readTab(const string& dataFN, const int label);
    static void writeTab(const string& dataFN, const SetHandler& norm,
                         const SetHandler& shuff);
    int const getLabel(int setPos);
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
