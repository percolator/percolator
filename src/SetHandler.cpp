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

#include "SetHandler.h"

SetHandler::SetHandler() {}

SetHandler::~SetHandler() {
  for (auto & subset : subsets) {
    if (subset != NULL) {
      delete subset;
    }
    subset = NULL;
  }
}

/**
 * Initialise without file input for @see PercolatorCInterface
 */
void SetHandler::filelessSetup(const unsigned int numFeatures,
                               const unsigned int numSpectra,
                               const set<int> labels) {
  for (auto label : labels) {
    DataSet* pSet = new DataSet();
    pSet->setLabel(label);
    pSet->setSize(numSpectra);
    pSet->initFeatureTables(numFeatures);
    subsets.push_back(pSet);
  }
}

/**
 * Insert DataSet object into this SetHandler
 * @param ds pointer to DataSet to be inserted
 */
void SetHandler::push_back_dataset( DataSet * ds ) {
  label2subset[ds->getLabel()] = subsets.size();  
  subsets.push_back(ds);
}

/**
 * Prints the results to a stream
 * @param test Scores object to be printed
 * @param myout stream to be printed to
 */
void SetHandler::print(Scores& test, int label, ostream& myout) {
  vector<ResultHolder> outList(0);
  for (auto & subset : subsets) {
    if (subset->getLabel() == label) {
      subset->print(test, outList);
    }
  }
  sort(outList.begin(), outList.end(), greater<ResultHolder> ());
  myout
      << "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
      << endl;
  for (const auto &psmResult : outList) {
    myout << psmResult << endl;
  }
}

/*const double* SetHandler::getFeatures(const int setPos, const int ixPos) const {
  return subsets[setPos]->getFeatures(ixPos);
}*/

int const SetHandler::getLabel(int setPos) {
  assert(setPos >= 0 && setPos < (signed int)subsets.size());
  return subsets[setPos]->getLabel();
}

void SetHandler::readTab(const string& dataFN) {
  if (VERB > 1) {
    cerr << "Reading Tab delimited input from datafile " << dataFN
        << endl;
  }
  
  // fill in the feature names from the first line
  ifstream dataStream(dataFN.c_str(), ios::out);
  if (!dataStream) {
    ostringstream temp;
    temp << "Error : Can not open file " << dataFN << endl;
    throw MyException(temp.str());
  }
  string tmp, line;
  getline(dataStream, line);
  istringstream iss(line);
  int skip = (DataSet::getCalcDoc() ? 2 : 0);
  while (iss.good()) {
    iss >> tmp;
    if (skip-- <= -2) { // removes enumerator, label and if present DOC features
      DataSet::getFeatureNames().insertFeature(tmp);
    }
  }
  
  // count number of features from first PSM
  getline(dataStream, line);
  unsigned int numFeatures = 0;
  iss.clear();
  iss.str(line);
  double a;
  iss >> tmp >> tmp; // remove id and label
  while (iss.good()) {
    iss >> a;
    ++numFeatures;
  }
  --numFeatures; // last one failed
  
  DataSet * targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(1);
  DataSet * decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(-1);
  
  if (!targetSet->initFeatures(numFeatures) || !decoySet->initFeatures(numFeatures)) {
    dataStream.close();
    throw MyException("Error : Reading tab file, too few features present.");
  }

  // read in the data
  string seq;
  int label;
  bool readError = false;
  dataStream.seekg(0, std::ios::beg);
  getline(dataStream, line); // skip over column names
  while (getline(dataStream, line)) {
    iss.clear();
    iss.str(line);
    iss >> tmp >> label;
    if (label == 1) {
      targetSet->readPsm(dataStream, line);
    } else if (label == -1) {
      decoySet->readPsm(dataStream, line);
    }
  }
  dataStream.close();

  
  push_back_dataset(targetSet);
  push_back_dataset(decoySet);
}

void SetHandler::writeTab(const string& dataFN) {
  ofstream dataStream(dataFN.data(), ios::out);
  dataStream << "SpecId\tLabel\t";
  if (DataSet::getCalcDoc()) {
    dataStream << "RT\tdM\t";
  }
  dataStream << DataSet::getFeatureNames().getFeatureNames(true)
      << "\tPeptide\tProteins" << endl;
  for (auto & subset : subsets) {
    subset->writeTabData(dataStream, subset->getLabel() == -1 ? "-1" : "1");
  }
  dataStream.close();
}
