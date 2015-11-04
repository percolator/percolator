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

SetHandler::SetHandler() : maxPSMs_(0u) {}

SetHandler::~SetHandler() {
  for (unsigned int ix = 0; ix < subsets_.size(); ix++) {
    if (subsets_[ix] != NULL) {
      delete subsets_[ix];
    }
    subsets_[ix] = NULL;
  }
}

/**
 * Gets the vector index of the DataSet matching the label
 * @param label DataSet label
 * @return index of matching DataSet
 */
unsigned int SetHandler::getSubsetIndexFromLabel(int label) {
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    if (subsets_[ix]->getLabel() == label) return ix;
  }
  ostringstream temp;
  temp << "Error: No DataSet found with label " << label << std::endl;
  throw MyException(temp.str());
}

/**
 * Insert DataSet object into this SetHandler
 * @param ds pointer to DataSet to be inserted
 */
void SetHandler::push_back_dataset( DataSet * ds ) {
  subsets_.push_back(ds);
}

/**
 * Prints the results to a stream
 * @param test Scores object to be printed
 * @param myout stream to be printed to
 */
void SetHandler::print(Scores& test, int label, ostream& myout) {
  vector<ResultHolder> outList(0);
  for (std::vector<DataSet*>::iterator it = subsets_.begin();
         it != subsets_.end(); ++it) {
    if ((*it)->getLabel() == label) {
      (*it)->print(test, outList);
    }
  }
  sort(outList.begin(), outList.end(), std::greater<ResultHolder> ());
  myout << "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n"
        << endl;
  for (std::vector<ResultHolder>::iterator it = outList.begin();
         it != outList.end(); ++it) {
    myout << *it << endl;
  }
}

void SetHandler::fillFeatures(vector<ScoreHolder> &scores, int label) {
  subsets_[getSubsetIndexFromLabel(label)]->fillFeatures(scores);
}

void SetHandler::normalizeFeatures(Normalizer*& pNorm) {
  std::vector<double*> featuresV, rtFeaturesV;
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    subsets_[ix]->fillFeatures(featuresV);
    subsets_[ix]->fillRtFeatures(rtFeaturesV);
  }
  pNorm = Normalizer::getNormalizer();

  pNorm->setSet(featuresV,
      rtFeaturesV,
      FeatureNames::getNumFeatures(),
      DataSet::getCalcDoc() ? RTModel::totalNumRTFeatures() : 0);
  pNorm->normalizeSet(featuresV, rtFeaturesV);
}

void SetHandler::setRetentionTime(map<int, double>& scan2rt) {
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    subsets_[ix]->setRetentionTime(scan2rt);
  }
}

/*const double* SetHandler::getFeatures(const int setPos, const int ixPos) const {
  return subsets_[setPos]->getFeatures(ixPos);
}*/

int const SetHandler::getLabel(int setPos) {
  assert(setPos >= 0 && setPos < (signed int)subsets_.size());
  return subsets_[setPos]->getLabel();
}

int SetHandler::readTab(istream& dataStream, SanityCheck*& pCheck) {
  if (!dataStream) {
    std::cerr << "ERROR: Cannot open data stream." << std::endl;
    return 0;
  }
  std::string psmLine, headerLine, defaultDirectionLine;
  
  getline(dataStream, headerLine); // line with feature names
  headerLine = rtrim(headerLine);
  if (headerLine.substr(0,5) == "<?xml") {
    std::cerr << "ERROR: Cannot read Tab delimited input from data stream.\n" << 
       "Input file seems to be in XML format, use the -k flag for XML input." << 
       std::endl;
    return 0;
  }
  
  // Checking for optional headers "ScanNr", "ExpMass" and "CalcMass"
  std::vector<OptionalField> optionalFields;
  int optionalFieldCount = getOptionalFields(headerLine, optionalFields);
  
  // parse second line for default direction
  getline(dataStream, defaultDirectionLine);
  defaultDirectionLine = rtrim(defaultDirectionLine);
  bool hasInitialValueRow = isDefaultDirectionLine(defaultDirectionLine);

  // count number of features from first PSM
  if (hasInitialValueRow) {
    getline(dataStream, psmLine);
  } else {
    psmLine = defaultDirectionLine;
  }
  psmLine = rtrim(psmLine);
  int numFeatures = getNumFeatures(psmLine, optionalFieldCount);
  
  // fill in the feature names from the header line
  FeatureNames& featureNames = DataSet::getFeatureNames();
  getFeatureNames(headerLine, numFeatures, optionalFieldCount, featureNames);
  
  // fill in the default weights if present
  std::vector<double> init_values;
  bool hasDefaultValues = false;
  if (hasInitialValueRow) {
    hasDefaultValues = getInitValues(defaultDirectionLine, optionalFieldCount, 
                                     init_values);
  }
  
  if (numFeatures < 1) {
    std::cerr << "ERROR: Reading tab file, too few features present." << std::endl;
    return 0;
  } else if (hasDefaultValues && init_values.size() > numFeatures) {
    std::cerr << "ERROR: Reading tab file, too many default values present." << std::endl;
    return 0;
  }

  // read in the data
  readPSMs(dataStream, psmLine, hasInitialValueRow, optionalFields);  
  //readSubsetPSMs(dataStream, psmLine, hasInitialValueRow, optionalFields, 1000u);
  
  pCheck = new SanityCheck();
  pCheck->checkAndSetDefaultDir();
  if (hasDefaultValues) pCheck->addDefaultWeights(init_values);
  return 1;
}

int SetHandler::getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields) {
  istringstream iss(headerLine);
  std::string tmp;
  iss >> tmp >> tmp; // discard id, label
  bool hasScannr = false;
  std::string optionalHeader;
  while (iss.good()) {
    iss >> optionalHeader;
    // transform to lower case for case insensitive matching
    std::transform(optionalHeader.begin(), optionalHeader.end(), 
                   optionalHeader.begin(), ::tolower);
    if (optionalHeader == "scannr") {
      optionalFields.push_back(SCANNR);
      hasScannr = true;
    } else if (optionalHeader == "expmass") {
      optionalFields.push_back(EXPMASS);
    } else if (optionalHeader == "calcmass") {
      optionalFields.push_back(CALCMASS);
    } else {
      break;
    }
  }
  if (!hasScannr) {
    cerr << "\nWARNING: Tab delimited input does not contain ScanNr column," <<
            "\n         scan numbers will be assigned automatically.\n" << endl;
  }
  return static_cast<int>(optionalFields.size());
}

bool SetHandler::isDefaultDirectionLine(const std::string& defaultDirectionLine) {
  istringstream iss(defaultDirectionLine);
  
  std::string psmid;
  iss >> psmid;

  std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
  return (psmid == "defaultdirection");
}

int SetHandler::getNumFeatures(const std::string& line, int optionalFieldCount) {
  istringstream iss(line);
  std::string tmp;
  iss >> tmp >> tmp; // remove id and label
  for (int i = 1; i <= optionalFieldCount; ++i) 
    iss >> tmp; // discard optional fields
  
  double a;
  int numFeatures = 0;
  iss >> a; // test third/fourth column
  while (iss.good()) {
    ++numFeatures;
    iss >> a;
  }
  
  if (DataSet::getCalcDoc()) numFeatures -= 2;
  return numFeatures;
}

void SetHandler::getFeatureNames(const std::string& headerLine, 
    int numFeatures, int optionalFieldCount, FeatureNames& featureNames) {
  istringstream iss(headerLine);
  int skip = 2 + optionalFieldCount + (DataSet::getCalcDoc() ? 2 : 0);
  int numFeatLeft = numFeatures;
  std::string tmp;
  while (iss.good()) {
    iss >> tmp;
    // removes enumerator, label and if present optional fields and DOC features
    if (skip-- <= 0 && numFeatLeft-- > 0) { 
      featureNames.insertFeature(tmp);
    }
  }
  
  featureNames.initFeatures(DataSet::getCalcDoc());
  assert(numFeatures == DataSet::getNumFeatures());
}

bool SetHandler::getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values) {
  istringstream iss(defaultDirectionLine);
  std::string tmp;
  iss >> tmp >> tmp; // remove id and label
  for (int i = 1; i <= optionalFieldCount + (DataSet::getCalcDoc() ? 2 : 0); ++i) 
    iss >> tmp; // discard optional fields
  
  bool hasDefaultValues = false;
  unsigned int ix = 0;
  double a;
  while (iss.good()) {
    iss >> a;
    if (a != 0.0) hasDefaultValues = true;
    if (VERB > 2) {
      std::cerr << "Initial direction for " << 
                   DataSet::getFeatureNames().getFeatureName(ix) << " is " << 
                   a << std::endl;
    }
    init_values.push_back(a);
    ix++;
  }
  return hasDefaultValues;
}

void SetHandler::writeTab(const string& dataFN, SanityCheck * pCheck) {
  ofstream dataStream(dataFN.data(), ios::out);
  dataStream << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\t";
  if (DataSet::getCalcDoc()) {
    dataStream << "RT\tdM\t";
  }
  dataStream << DataSet::getFeatureNames().getFeatureNames(true)
      << "\tPeptide\tProteins" << std::endl;
  vector<double> initial_values = pCheck->getDefaultWeights();
  if (initial_values.size() > 0) {
    dataStream << "DefaultDirection\t-\t-\t-\t-";
    if (DataSet::getCalcDoc()) {
      dataStream << "\t-\t-";
    }
    for (size_t i = 0; i < initial_values.size(); ++i) {
      dataStream << '\t' << initial_values[i];
    }
    dataStream << std::endl;
  }
  for (std::vector<DataSet*>::iterator it = subsets_.begin();
         it != subsets_.end(); ++it) {
    (*it)->writeTabData(dataStream);
  }
}

void SetHandler::readPSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields) {
  DataSet* targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(1);
  DataSet* decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(-1);
  
  if (maxPSMs_ > 0u) {
    std::priority_queue<PSMDescriptionPriority> subsetPSMs;
    std::map<ScanId, size_t> scanIdLookUp;
    unsigned int lineNr = (hasInitialValueRow ? 3u : 2u);
    do {
      psmLine = rtrim(psmLine);
      ScanId scanId = getScanId(psmLine, optionalFields, lineNr);
      size_t randIdx;
      if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
        randIdx = scanIdLookUp[scanId];
      } else {
        randIdx = Scores::lcg_rand();
        scanIdLookUp[scanId] = randIdx;
      }
      
      unsigned int upperLimit = UINT_MAX;
      if (subsetPSMs.size() < maxPSMs_ || randIdx < upperLimit) {
        PSMDescriptionPriority psmPriority;
        psmPriority.psm = new PSMDescription();
        bool readProteins = false;
        psmPriority.label = DataSet::readPsm(psmLine, lineNr, optionalFields, 
                                             readProteins, psmPriority.psm);
        psmPriority.priority = randIdx;
        subsetPSMs.push(psmPriority);
        if (subsetPSMs.size() > maxPSMs_) {
          PSMDescriptionPriority del = subsetPSMs.top();
          upperLimit = del.priority;
          deletePSMPointer(del.psm);
          subsetPSMs.pop();
        }
      }
      ++lineNr;
    } while (getline(dataStream, psmLine));
    
    while (!subsetPSMs.empty()) {
      PSMDescriptionPriority psmPriority = subsetPSMs.top();
      if (psmPriority.label == 1) {
        targetSet->registerPsm(psmPriority.psm);
      } else if (psmPriority.label == -1) {
        decoySet->registerPsm(psmPriority.psm);
      } else {
        std::cerr << "Warning: the PSM " << psmPriority.psm->id
            << " has a label not in {1,-1} and will be ignored." << std::endl;
        deletePSMPointer(psmPriority.psm);
      }
      subsetPSMs.pop();
    }
  } else {
    unsigned int targetIdx = 0u, decoyIdx = 0u, lineNr = (hasInitialValueRow ? 3u : 2u);
    do {
      psmLine = rtrim(psmLine);
      int label = getLabel(psmLine, lineNr);
      if (label == 1) {
        targetSet->readPsm(psmLine, lineNr, optionalFields);
      } else if (label == -1) {
        decoySet->readPsm(psmLine, lineNr, optionalFields);
      } else {
        std::cerr << "Warning: the PSM on line " << lineNr
            << " has a label not in {1,-1} and will be ignored." << std::endl;
      }
      ++lineNr;
    } while (getline(dataStream, psmLine));
  }
  
  push_back_dataset(targetSet);
  push_back_dataset(decoySet);
}

int SetHandler::getLabel(const std::string& psmLine, unsigned int lineNr) {
  istringstream iss(psmLine);
  std::string psmid;
  int label;
  iss >> psmid >> label;
  if (!iss.good()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read PSMid or label." << std::endl;
    throw MyException(temp.str());
  }
  return label;
}

ScanId SetHandler::getScanId(const std::string& psmLine, 
    std::vector<OptionalField>& optionalFields, unsigned int lineNr) {
  ScanId scanId;
  
  std::istringstream buff(psmLine);
  std::string tmp;
  int label;
  buff >> tmp >> label; // read PSMid and get rid of label
  if (!buff.good()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read PSMid or label." << std::endl;
    throw MyException(temp.str());
  }
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        buff >> scanId.first;
        if (!buff.good()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number on line " 
              << lineNr << ". Check if scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        buff >> scanId.second;
        if (!buff.good()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading experimental mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case CALCMASS: {
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) scanId.first = lineNr;
  
  return scanId;
}

void SetHandler::deletePSMPointer(PSMDescription* psm) {
  if (psm->features) {
    delete[] psm->features;
    psm->features = NULL;
  }
  if (psm->retentionFeatures) {
    delete[] psm->retentionFeatures;
    psm->retentionFeatures = NULL;
  }
  delete psm;
}

std::string& SetHandler::rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}
