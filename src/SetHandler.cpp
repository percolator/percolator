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

SetHandler::SetHandler() {
  n_examples = 0;
  labels = NULL;
  c_vec = NULL;
}

SetHandler::~SetHandler() {
  if (labels) {
    delete[] labels;
  }
  labels = NULL;
  if (c_vec) {
    delete[] c_vec;
  }
  c_vec = NULL;
  for (unsigned int ix = 0; ix < subsets.size(); ix++) {
    if (subsets[ix] != NULL) {
      delete subsets[ix];
    }
    subsets[ix] = NULL;
  }
}

void SetHandler::filelessSetup(const unsigned int numFeatures,
                               const unsigned int numSpectra,
                               const int label) {
  DataSet* pSet = new DataSet();
  pSet->setLabel(label);
  pSet->initFeatureTables(numFeatures);
  subsets.push_back(pSet);
  n_examples = numSpectra;
}




void SetHandler::push_back_dataset( DataSet * ds ) {
    subsets.push_back(ds);
}


void SetHandler::print(Scores& test, ostream& myout) {
  vector<ResultHolder> outList(0);
  for (unsigned int setPos = 0; setPos < subsets.size(); setPos++) {
    subsets[setPos]->print(test, outList);
  }
  sort(outList.begin(), outList.end(), greater<ResultHolder> ());
  vector<ResultHolder>::iterator it = outList.begin();
  myout
      << "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
      << endl;
  for (; it != outList.end(); ++it) {
    myout << *it << endl;
  }
}

void SetHandler::generateTrainingSet(const double fdr, const double cpos,
                                     const double cneg, Scores& sc) {
  double tp = 0, fp = 0;
  unsigned int ix = 0;
  examples.clear();
  bool underCutOff = true;
  vector<ScoreHolder>::const_iterator it;
  for (it = sc.begin(); it != sc.end(); it++) {
    if (it->label == -1) {
      fp++;
    } else {
      tp++;
    }
    if (underCutOff && fdr < (fp / (tp + fp))) {
      underCutOff = false;
    }
    if (it->label == -1 || underCutOff) {
      examples.push_back(it->pPSM->features);
      labels[ix] = it->label;
      c_vec[ix++] = (it->label != -1 ? cpos : cneg);
    }
  }
}

PSMDescription* SetHandler::getNext(int& setPos, int& ixPos) {
  PSMDescription* features = subsets[setPos]->getNext(ixPos);
  if (features) {
    return features;
  }
  if (++setPos >= ((signed int)subsets.size())) {
    return NULL;
  }
  ixPos = -1;
  return subsets[setPos]->getNext(ixPos);
}

/*const double* SetHandler::getFeatures(const int setPos, const int ixPos) const {
  return subsets[setPos]->getFeatures(ixPos);
}*/

int const SetHandler::getLabel(int setPos) {
  assert(setPos >= 0 && setPos < (signed int)subsets.size());
  return subsets[setPos]->getLabel();
}

void SetHandler::setSet() {
  n_examples = 0;
  int i = 0, j = -1;
  while (getNext(i, j)) {
    n_examples++;
  }
  if (!labels) {
    labels = new double[n_examples];
  }
  if (!c_vec) {
    c_vec = new double[n_examples];
  }
  if (VERB > 3) {
    cerr << "Set up a SetHandler with " << subsets.size()
        << " DataSet:s and " << n_examples << " examples" << endl;
    if (VERB > 4) {
      for (unsigned int i = 0; i < subsets.size(); i++) {
        cerr << "First 10 lines of " << i + 1 << " set with "
            << subsets[i]->getLabel() << " label" << endl;
        subsets[i]->print_10features();
      }
    }
  }
}

void SetHandler::readTab(const string& dataFN, const int setLabel) {
  if (VERB > 1) {
    cerr << "Reading Tab delimetered input from datafile " << dataFN
        << endl;
  }
  ifstream labelStream(dataFN.c_str(), ios::out);
  if (!labelStream) {
    ostringstream temp;
    temp << "Error : Can not open file " << dataFN << endl;
    throw MyException(temp.str());
  }
  vector<unsigned int> ixs;
  ixs.clear();
  string tmp, line;
  int label;
  unsigned int ix = 0;
  getline(labelStream, tmp); // Id row
  while (true) {
    labelStream >> tmp >> label;
    getline(labelStream, tmp); // read rest of line
    if (!labelStream) {
      break;
    }
    if (label == setLabel) {
      ixs.push_back(ix);
    }
    ++ix;
  }
  labelStream.close();
  ifstream dataStream(dataFN.c_str(), ios::out);
  if (!dataStream) {
    ostringstream temp;
    temp << "Error : Can not open file " << dataFN << endl;
    throw MyException(temp.str());
  }
  dataStream >> tmp >> tmp;
  dataStream.get(); // removed enumrator, label and tab
  getline(dataStream, line);
  istringstream iss(line);
  int skip = (DataSet::getCalcDoc() ? 2 : 0);
  while (iss.good()) {
    iss >> tmp;
    if (skip-- <= 0) {
      DataSet::getFeatureNames().insertFeature(tmp);
    }
  }
  DataSet* theSet = new DataSet();
  theSet->setLabel(setLabel > 0 ? 1 : -1);
  cerr << "Point readTabData" << endl;
  theSet->readTabData(dataStream, ixs);
  cerr << "Point dataStream.close()" << endl;
  dataStream.close();
  subsets.push_back(theSet);
  setSet();
}

void SetHandler::writeTab(const string& dataFN, const SetHandler& norm,
                          const SetHandler& shuff) {
  ofstream dataStream(dataFN.data(), ios::out);
  dataStream << "SpecId\tLabel\t";
  if (DataSet::getCalcDoc()) {
    dataStream << "RT\tdM\t";
  }
  dataStream << DataSet::getFeatureNames().getFeatureNames(true)
      << "\tPeptide\tProteins" << endl;
  string str;
  for (int setPos = 0; setPos < (signed int)norm.subsets.size(); setPos++) {
    norm.subsets[setPos]->writeTabData(dataStream,
                                       norm.subsets[setPos]->getLabel()
                                           == -1 ? "-1" : "1");
  }
  for (int setPos = 0; setPos < (signed int)shuff.subsets.size(); setPos++) {
    shuff.subsets[setPos]->writeTabData(dataStream,
                                        shuff.subsets[setPos]->getLabel()
                                            == -1 ? "-1" : "1");
  }
  dataStream.close();
}
