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
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    if (psm->hasSpectrumFileName()) {
      out << '\t' << psm->getSpectrumFileName();
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << psm->getFullPeptide() << '\t';
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
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
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

/**
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
void DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, FeatureMemoryPool& featurePool,
    std::string decoyPrefix) { 
  PSMDescription* myPsm = NULL;
  bool readProteins = true;
  readPsm(line, lineNr, optionalFields, readProteins, myPsm, featurePool, decoyPrefix);
  registerPsm(myPsm);
}

bool TabFileValidator::decoyWarningTripped = false;

int DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, std::string decoyPrefix) {
  TabReader reader(line);
  std::string tmp;
  
  myPsm = new PSMDescription();
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
          ostringstream temp;
          temp << "ERROR: The parser tries to allocate a deprecated delta mass field for the PSM " 
              << myPsm->getId() << std::endl;
          throw MyException(temp.str());
        break;
      } case FILENAME: {
        myPsm->setSpectrumFileName(reader.readString());
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
  myPsm->setPeptide(peptide_seq);
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM " << myPsm->getId() 
      << ". Check if a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (ProteinProbEstimator::getCalcProteinLevelProb()) {
    // MT: we only need the peptide sequences to be well formatted if protein inference is applied
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
    if (PSMDescription::getProteinNameSeparator() == "\t") {
      while (!reader.error()) {
        std::string tmp = reader.readString();
        if (tmp.size() > 0) proteins.push_back(tmp);
      }
    } else {
      std::string names = reader.readString();
      size_t pos = 0;
      std::string token;
      while ((pos = names.find(PSMDescription::getProteinNameSeparator())) != std::string::npos) {
        token = names.substr(0, pos);
        if (token.size() > 0) proteins.push_back(token);
        names.erase(0, pos + PSMDescription::getProteinNameSeparator().length());
      }
      if (names.size() > 0) proteins.push_back(names);
    }
    proteins.swap(myPsm->proteinIds); // shrink to fit
  }

  if (label == -1) {
    for(auto const& value: myPsm->proteinIds) { 
      std::string token = value.substr(0, value.find("_") + 1);
      if (token != decoyPrefix && VERB > 1 && TabFileValidator::decoyWarningTripped == false) {
        std::cerr << "Warning: Set decoy prefix don't match" << std::endl;
        TabFileValidator::decoyWarningTripped = true;
      }
    }
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
  psms_.push_back(myPsm);
}

bool TabFileValidator::isTabFile(std::string file_name) {
  // open C++ stream to file
  std::ifstream file(file_name.c_str());
  // file not opened, return false
  if(!file.is_open()) return false;
  // read a line from the file       
  std::string wtf;
  std::istream &in= std::getline(file, wtf);
  // unable to read the line, return false
  if(!in) {
    if (VERB > 0) {
      std::cerr << "Cannot read" << file_name << "!" << std::endl;
    }
    return false;
  }
  // try to find a '\t', return true if '\t' is found within the string
  bool isTab = std::find(wtf.begin(), wtf.end(), '\t')!= wtf.end();
  if (!isTab && VERB > 0) {
    std::cerr << file_name << " is not comma delimited!\n" << std::endl;
  }
  return isTab;
}

bool TabFileValidator::isTabFiles(std::vector<std::string> files) {
  for (const string &file : files) {
    if(!isTabFile(file)) {
      return false;
    } 
  };
  return true;
}

std::string TabFileValidator::getDecoyPrefix(std::vector<std::string> fileList) {
  std::string decoy_prefix;
  for(const string &file : fileList) {
      decoy_prefix = detectDecoyPrefix(file);
      break;
  };
  return decoy_prefix;
}

void TabFileValidator::getProteinIndex(std::string file_name, int* proteinIndex, int* labelIndex) {
  // open C++ stream to file
  std::ifstream file(file_name.c_str());

  // file not opened, return false
  if(!file.is_open()) {
    return ;
  }

  int columnIndex = 0;
  std::string headerRow;
  getline(file, headerRow);
  TabReader reader(headerRow);
    while (!reader.error()) {
      std::string optionalHeader = reader.readString();
      std::transform(optionalHeader.begin(), optionalHeader.end(), optionalHeader.begin(),
    [](unsigned char c){ return std::tolower(c); });
      if (optionalHeader.find("proteins") != std::string::npos) { 
        *proteinIndex = columnIndex;
      } else if (optionalHeader.find("label") != std::string::npos){
        *labelIndex = columnIndex;
      }
      columnIndex++;
    }
}

std::string TabFileValidator::LongestCommonSubsequence(vector<string> arr)
{
    string res = "";
    // Determine size of the array
    int n = arr.size();
    if (n<=0) { return res; } // Check for limited size
 
    // Take first word from array as reference
    string s = arr[0];
    if (n==1) { return s; } // Check for limited size
    int len = s.length();
  
    for (int i = 0; i < len; i++) {
        for (int j = i + 1; j <= len; j++) {
            // generating all possible substrings
            // of our reference string arr[0] i.e s
            string stem = s.substr(i, j);
            int k = 1;
            for (; k < n; k++) {
                // Check if the generated stem is
                // common to all words
                if (arr[k].find(stem) == std::string::npos)
                    break;
            }
 
            // If current substring is present in
            // all strings and its length is greater
            // than current result
            if (k == n && res.length() < stem.length())
                res = stem;
        }
    }
 
    return res;
}

std::string TabFileValidator::findDecoyPrefix(std::string file_name, int proteinIndex, int labelIndex) {
  // open C++ stream to file
  std::ifstream file(file_name.c_str());

  // file not opened, return false
  if(!file.is_open()) return "error";
  // read a line from the file  

  /* Skip header */
  std::string nextRow;
  
  getline(file, nextRow);

  std::vector<std::string> proteinNames;
  while(getline(file, nextRow)) {
    TabReader readerRow(nextRow);
    int col = 0;
    bool isDecoy = false;
    while (!readerRow.error()) {
      std::string value = readerRow.readString();
      if (col == labelIndex && value == "-1") {
        isDecoy = true;
      }
      if (col == proteinIndex && isDecoy) {
        proteinNames.push_back(value);
        
      }
      col++;
    }
  }

  /* Shuffle protein ids  */
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(proteinNames), std::end(proteinNames), rng);
  /* Check prefix for 10% of all protein ids  */
  int n_elements = round(0.1 * proteinNames.size());
  proteinNames.resize(n_elements);
  
  std::string prefix = LongestCommonSubsequence(proteinNames);
  size_t loc = prefix.find("_");
  if (loc != std::string::npos) {
    prefix = prefix.substr(0, loc + 1);
  }

  if (VERB > 0) {
    std::cerr << "Protein decoy-prefix used is " << prefix << std::endl;  
  }
  
  return prefix;
}



std::string TabFileValidator::detectDecoyPrefix(std::string file_name) {
  int proteinIndex=-1;
  int labelIndex=-1;

  getProteinIndex(file_name, &proteinIndex, &labelIndex);
  if (proteinIndex==-1 || labelIndex==-1) {
    if (VERB > 0) {
        std::cerr <<  proteinIndex << " " << labelIndex << std::endl;
 
        std::cerr << "Couldn't find 'Proteins' column in tab-file" << std::endl;
    }
    return "error";
  }
  return findDecoyPrefix(file_name, proteinIndex, labelIndex);
}

bool TabFileValidator::validateTabFiles(std::vector<std::string> files, std::string* decoy_prefix) {
  std::string tmpDecoyPrefix = getDecoyPrefix(files);
  *decoy_prefix = tmpDecoyPrefix;

  if (tmpDecoyPrefix=="error") {
    return false;
  }
  if (!isTabFiles(files)) {
    return false;
  }
  return true;
}
