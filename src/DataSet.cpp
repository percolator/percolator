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

bool DataSet::isotopeMass = false;
bool DataSet::calcDOC = false;
const string DataSet::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";
string DataSet::ptmAlphabet = "#*@";
FeatureNames DataSet::featureNames;

DataSet::DataSet() : numSpectra(0) {}

DataSet::~DataSet() {
  BOOST_FOREACH (PSMDescription * psm, psms)
  {
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
  if (calcDOC) {
    nf -= DescriptionOfCorrect::numDOCFeatures();
  }
  BOOST_FOREACH (PSMDescription * pPSM, psms) {
    double* featureRow = pPSM->features;
    out << pPSM->id << '\t' << label << '\t' << pPSM->scan;
    if (calcDOC) {
      out << '\t' << pPSM->getUnnormalizedRetentionTime() << '\t'
          << pPSM->massDiff;
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << pPSM->peptide;
    BOOST_FOREACH (const std::string &proteinId, pPSM->proteinIds) {
      out << '\t' << proteinId;
    }
    out << endl;
  }
  return true;
}

void DataSet::print_features() {
  for (int i = 0; i < getSize(); i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << psms[i]->features[j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features() {
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (int i = 0; i < 10; i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

void DataSet::print(Scores& test, vector<ResultHolder> &outList) {
  ostringstream out;
  BOOST_FOREACH (PSMDescription * psm, psms) {
    ScoreHolder* pSH = test.getScoreHolder(psm->features);
    if (pSH == NULL) {
      continue;
    }
    BOOST_FOREACH (const std::string &proteinId, psm->proteinIds) {
      out << "\t" << proteinId;
    }
    outList.push_back(ResultHolder(pSH->score, psm->q, psm->pep, psm->id, psm->peptide, out.str()));
    out.str("");
  }
}

// TODO: find a way to make these four functions generic
void DataSet::fillFeatures(vector<ScoreHolder> &scores) {
  BOOST_FOREACH (PSMDescription * psm, psms)
    scores.push_back(ScoreHolder(.0, label, psm));
}
    
void DataSet::fillFeaturesPeptide(vector<ScoreHolder> &scores) {
  BOOST_FOREACH (PSMDescription * psm, psms)
    scores.push_back(ScoreHolderPeptide(.0, label, psm));
}

void DataSet::fillFeatures(vector<double*> &features) {
  BOOST_FOREACH (PSMDescription * psm, psms)
    features.push_back(psm->features);
}

void DataSet::fillRtFeatures(vector<double*> &rtFeatures) {
  double* features;
  BOOST_FOREACH (PSMDescription * psm, psms)
    if ((features = psm->retentionFeatures))
        rtFeatures.push_back(features);
}

/*
double DataSet::isPngasef(const string& peptide) {
  bool isDecoy;
  switch (label) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("ERROR : class DataSet has not been initiated\
		to neither target nor decoy label\n");}
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
void DataSet::readPsm(ifstream & dataStream, const std::string line) {
  istringstream buff(line);
  string tmp;
  unsigned int numFeatures = FeatureNames::getNumFeatures();
  
  PSMDescription *myPsm = new PSMDescription();
  buff >> myPsm->id;
  buff >> tmp; // get rid of label
  buff >> myPsm->scan;
  if (calcDOC) {
    numFeatures -= DescriptionOfCorrect::numDOCFeatures();
    buff >> myPsm->retentionTime;
    buff >> myPsm->massDiff;
  }
  double *featureRow = new double[numFeatures];
  myPsm->features = featureRow;
  for (register unsigned int j = 0; j < numFeatures; j++) {
    buff >> featureRow[j];
  }  
  std::string peptide_seq = "";
  buff >> peptide_seq;
  myPsm->peptide = peptide_seq;
  
  // do some error checking
  if (!buff.good()) {
    dataStream.close();
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM with id " << myPsm->id << ". Check if\
    the line is formatted correctly." << std::endl;
    throw MyException(temp.str());
  } else if (peptide_seq.size() < 5) {
    dataStream.close();
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq << "\
    with PSM id " << myPsm->id << " is too short." << std::endl;
    throw MyException(temp.str());
  } else if (peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.') {
    dataStream.close();
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq << "\
    with PSM id " << myPsm->id << " does not contain one or two of its flanking amino acids." << std::endl;
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
    if (aaAlphabet.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

unsigned int DataSet::cntPTMs(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (ptmAlphabet.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

void DataSet::registerPsm(PSMDescription * myPsm) {
  bool isDecoy;
  switch (label) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("Error : Reading PSM, class DataSet has not been initiated\
		to neither target nor decoy label\n");}
  }
  
  if (calcDOC) {
    int featureNum = featureNames.getDocFeatNum();
    myPsm->retentionFeatures = new double[RTModel::totalNumRTFeatures()];
    DescriptionOfCorrect::calcRegressionFeature(*myPsm);
    myPsm->features[featureNum++] = abs( myPsm->pI - 6.5);
    myPsm->features[featureNum++] = abs( myPsm->massDiff);
    myPsm->features[featureNum++] = 0;
    myPsm->features[featureNum++] = 0;
  }
  psms.push_back(myPsm);
  ++numSpectra;
}
