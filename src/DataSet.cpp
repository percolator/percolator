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

#ifdef WIN32
#include <float.h>
#define isfinite _finite
#else
#include <cmath>
#endif

bool DataSet::calcDOC_ = false;
FeatureNames DataSet::featureNames_;

DataSet::DataSet() {}

DataSet::~DataSet() {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription::deletePtr(*it);
  }
}

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  if (calcDOC_) {
    nf -= static_cast<unsigned int>(DescriptionOfCorrect::numDOCFeatures());
  }
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    if (calcDOC_) {
      out << '\t' << psm->getUnnormalizedRetentionTime() << '\t'
          << psm->getMassDiff();
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
  for (std::size_t i = 0; i < getSize(); i++) {
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << psms_[i]->features[j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features() {
  
  for (std::size_t i = 0; i < 10; i++) {
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms_[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

// TODO: find a way to make these three functions generic
void DataSet::fillScores(std::vector<ScoreHolder>& scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm));
  }
}

void DataSet::fillFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

void DataSet::fillDOCFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features + featureNames_.getDocFeatNum());
  }
}

void DataSet::fillRtFeatures(std::vector<double*>& rtFeatures) {
  double* features;
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    if ((features = psm->getRetentionFeatures()))
      rtFeatures.push_back(features);
  }
}

/**
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
void DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, FeatureMemoryPool& featurePool) { 
  PSMDescription* myPsm = NULL;
  bool readProteins = true;
  readPsm(line, lineNr, optionalFields, readProteins, myPsm, featurePool);
  registerPsm(myPsm);
}

int DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool) {
  TabReader reader(line);
  std::string tmp;
  
  if (calcDOC_) {
    myPsm = new PSMDescriptionDOC();
  } else {
    /* myPsm = new PSMDescription(); */
    myPsm = new PSMDescriptionDOC();
  }
  myPsm->setId(reader.readString());
  int label = reader.readInt();
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        myPsm->scan = static_cast<unsigned int>(reader.readInt());
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number of PSM " 
              << myPsm->getId() << ". Check if scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        myPsm->expMass = reader.readDouble();
        break;
      } case CALCMASS: {
        myPsm->calcMass = reader.readDouble();
        break;
      } case RETTIME: {
        myPsm->setRetentionTime(reader.readDouble());
        break;
      } case DELTAMASS: {
        myPsm->setMassDiff(reader.readDouble());
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
  
  unsigned int numFeatures = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  if (calcDOC_) {
    numFeatures -= static_cast<unsigned int>(DescriptionOfCorrect::numDOCFeatures());
  }
  double* featureRow = featurePool.allocate();
  myPsm->features = featureRow;
  for (register unsigned int j = 0; j < numFeatures; j++) {
    featureRow[j] = reader.readDouble();
    if (!isfinite(featureRow[j])) {
      ostringstream oss;
      oss << "ERROR: Reached strange feature with val=" << featureRow[j]
           << " col=" << j << " for PSM with id " << myPsm->getId() << endl;
      if (NO_TERMINATE) {
        cerr << oss.str();
        std::cerr << "No-terminate flag set: setting value to 0 and ignoring the error." << std::endl;
        featureRow[j] = 0.0;
      } else {
        throw MyException(oss.str() + "Terminating.\n");
      }
    }
  }
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in feature vector of PSM " 
      << myPsm->getId() << ". Check if there are enough features on this line and "
      << "if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  std::string peptide_seq = reader.readString();
  myPsm->peptide = peptide_seq;
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM " << myPsm->getId() 
      << ". Check if a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (calcDOC_ || ProteinProbEstimator::getCalcProteinLevelProb()) {
    // MT: we only need the peptide sequences to be well formatted if DOC features 
    // are calculated, or if protein inference is applied
    if (peptide_seq.size() < 5) {
      ostringstream temp;
      temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq 
        << " with PSM id " << myPsm->getId() << " is too short." << std::endl;
      throw MyException(temp.str());
    } else if (peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.') {
      ostringstream temp;
      temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq 
        << " with PSM id " << myPsm->getId() << " does not contain one or two of its"
        << " flanking amino acids." << std::endl;
      throw MyException(temp.str());
    }
  }
  
  if (readProteins) {
    std::vector<std::string> proteins;
    while (!reader.error()) {
      std::string tmp = reader.readString();
      if (tmp.size() > 0) proteins.push_back(tmp);
    }
    proteins.swap(myPsm->proteinIds); // shrink to fit
  }
  
  return label;
}

void DataSet::registerPsm(PSMDescription* myPsm) {
  switch (label_) {
    case 1: { break; };
    case -1: { break; };
    default:  { throw MyException("ERROR : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label\n");}
  }
  
  if (calcDOC_) {
    myPsm->setRetentionFeatures(new double[RTModel::totalNumRTFeatures()]());
    DescriptionOfCorrect::calcRegressionFeature(myPsm);
  }
  psms_.push_back(myPsm);
}
