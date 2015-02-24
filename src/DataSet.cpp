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

#include "DataSet.h"

bool DataSet::calcDOC_ = false;
const string DataSet::aaAlphabet_ = "ACDEFGHIKLMNPQRSTVWY";
string DataSet::ptmAlphabet_ = "#*@";
FeatureNames DataSet::featureNames_;

DataSet::DataSet() : numSpectra_(0) {}

DataSet::~DataSet() {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    if(psm->features) {
      delete[] psm->features;
      psm->features = NULL;
    }
    if(psm->retentionFeatures) {
      delete[] psm->retentionFeatures;
      psm->retentionFeatures = NULL;
    }
    delete psm;
    psm = NULL;
  }
}

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = FeatureNames::getNumFeatures();
  if (calcDOC_) {
    nf -= DescriptionOfCorrect::numDOCFeatures();
  }
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->id << '\t' << label_ << '\t' << psm->scan;
    if (calcDOC_) {
      out << '\t' << psm->getUnnormalizedRetentionTime() << '\t'
          << psm->massDiff;
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << psm->peptide;
    psm->printProteins(out);
    out << endl;
  }
  return true;
}

void DataSet::print_features() {
  for (int i = 0; i < getSize(); i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << psms_[i]->features[j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features() {
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (int i = 0; i < 10; i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms_[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

void DataSet::print(Scores& test, vector<ResultHolder> &outList) {
  ostringstream out;
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    ScoreHolder* pSH = test.getScoreHolder(psm->features);
    if (pSH == NULL) {
      continue;
    }
    psm->printProteins(out);
    outList.push_back(ResultHolder(pSH->score, pSH->q, pSH->pep, psm->id, psm->peptide, out.str()));
    out.str("");
  }
  test.resetScoreMap();
}

// TODO: find a way to make these three functions generic
void DataSet::fillFeatures(vector<ScoreHolder> &scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm));
  }
}

void DataSet::fillFeatures(vector<double*> &features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

void DataSet::fillRtFeatures(vector<double*> &rtFeatures) {
  double* features;
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    if ((features = psm->retentionFeatures))
        rtFeatures.push_back(features);
  }
}

/*
double DataSet::isPngasef(const string& peptide) {
  bool isDecoy;
  switch (label_) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("ERROR : class DataSet has not been initiated\
    to neither target nor decoy label_\n");}
  }
  return isPngasef( peptide, isDecoy);
}
//
double DataSet::isPngasef(const string& peptide, bool isDecoy ) {
  size_t next_pos = 0, pos;
  while ((pos = peptide.find("N*", next_pos)) != string::npos) {
    next_pos = pos + 1;
    if (! isDecoy) {
      pos += 3;
      if (peptide[pos] == '#') {
        pos += 1;
      }
    } else {
      pos -= 2;
      if (peptide[pos] == '#') {
        pos -= 1;
      }
    }
    if (peptide[pos] == 'T' || peptide[pos] == 'S') {
      return 1.0;
    }
  }
  return 0.0;
}
*/

/**
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
void DataSet::readPsm(const std::string line, const unsigned int lineNr, 
                      const std::vector<OptionalField>& optionalFields) {
  istringstream buff(line);
  string tmp;
  unsigned int numFeatures = FeatureNames::getNumFeatures();
  
  PSMDescription *myPsm = new PSMDescription();
  buff >> myPsm->id >> tmp; // read PSMid and get rid of label_
  if (!buff.good()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr << \
        ". Could not read PSMid or label_." << std::endl;
    throw MyException(temp.str());
  }
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        buff >> myPsm->scan;
        if (!buff.good()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number of PSM with id " << myPsm->id << \
              ". Check if the scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        buff >> myPsm->expMass;
        break;
      } case CALCMASS: {
        buff >> myPsm->calcMass;
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) myPsm->scan = lineNr;
  
  if (calcDOC_) {
    numFeatures -= DescriptionOfCorrect::numDOCFeatures();
    buff >> myPsm->retentionTime;
    buff >> myPsm->massDiff;
  }
  double *featureRow = new double[FeatureNames::getNumFeatures()];
  myPsm->features = featureRow;
  for (register unsigned int j = 0; j < numFeatures; j++) {
    buff >> featureRow[j];
  }
  if (!buff.good()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in the feature vector of PSM with id " << myPsm->id << \
     ". Check if there are enough features on this line and if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  std::string peptide_seq = "";
  buff >> peptide_seq;
  myPsm->peptide = peptide_seq;
  if (!buff.good()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM with id " << myPsm->id << ". Check if" << \
      " a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (peptide_seq.size() < 5) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq << \
      " with PSM id " << myPsm->id << " is too short." << std::endl;
    throw MyException(temp.str());
  } else if (peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.') {
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq << \
      " with PSM id " << myPsm->id << " does not contain one or two of its flanking amino acids." << std::endl;
    throw MyException(temp.str());
  }
  
  while (!!buff) {
    buff >> tmp;
    if (tmp.size() > 0) {
      myPsm->proteinIds.insert(tmp);
    }
  }
  
  registerPsm(myPsm);
}

unsigned int DataSet::peptideLength(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (aaAlphabet_.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

unsigned int DataSet::cntPTMs(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (ptmAlphabet_.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

void DataSet::registerPsm(PSMDescription * myPsm) {
  bool isDecoy;
  switch (label_) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("Error : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label_\n");}
  }
  
  if (calcDOC_) {
    int featureNum = featureNames_.getDocFeatNum();
    myPsm->retentionFeatures = new double[RTModel::totalNumRTFeatures()];
    DescriptionOfCorrect::calcRegressionFeature(*myPsm);
    myPsm->features[featureNum++] = abs( myPsm->pI - 6.5);
    myPsm->features[featureNum++] = abs( myPsm->massDiff);
    myPsm->features[featureNum++] = 0;
    myPsm->features[featureNum++] = 0;
  }
  psms_.push_back(myPsm);
  ++numSpectra_;
}
