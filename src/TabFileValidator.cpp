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

#include "TabFileValidator.h"

bool TabFileValidator::isTabFile(std::string fileName) {
  std::ifstream file(fileName.c_str());
  
  if (!file.is_open()) return false;
  
  std::string tmp;
  std::istream &in = std::getline(file, tmp);
  
  if (!in) {
    if (VERB > 0) {
      std::cerr << "Cannot read " << fileName << "!" << std::endl;
    }
    return false;
  }
  // try to find a '\t', return true if '\t' is found within the string
  bool isTab = std::find(tmp.begin(), tmp.end(), '\t')!= tmp.end();
  if (!isTab && VERB > 0) {
    std::cerr << fileName << " is not tab delimited!\n" << std::endl;
  }
  return isTab;
}

bool TabFileValidator::isTabFiles(std::vector<std::string> fileNames) {
  for (const string &fileName : fileNames) {
    if (!isTabFile(fileName)) {
      return false;
    } 
  };
  return true;
}

std::string TabFileValidator::getDecoyPrefix(std::vector<std::string> fileNames) {
  std::string decoyPrefix;
  for (const string &fileName : fileNames) {
      decoyPrefix = detectDecoyPrefix(fileName);
      break;
  };
  return decoyPrefix;
}

void TabFileValidator::getProteinAndLabelColumnIndices(std::string fileName, int &proteinIndex, int &labelIndex) {
  // open C++ stream to file
  std::ifstream file(fileName.c_str());

  // file not opened, return false
  if(!file.is_open()) {
    return ;
  }

  std::string headerRow;
  getline(file, headerRow);

  TabReader reader(headerRow);
  int columnIndex = 0;
  while (!reader.error()) {
    std::string optionalHeader = reader.readString();
    std::transform(
      optionalHeader.begin(),
      optionalHeader.end(),
      optionalHeader.begin(),
      [](unsigned char c){ return std::tolower(c); }
    );
    if (optionalHeader.find("proteins") != std::string::npos) { 
      proteinIndex = columnIndex;
    } else if (optionalHeader.find("label") != std::string::npos){
      labelIndex = columnIndex;
    }
    columnIndex++;
  }
}

std::string TabFileValidator::getLongestCommonPrefix(std::vector<std::string> strings) {
  size_t numStrings = strings.size();
  if (numStrings == 0) {
    return "";
  }

  std::string referenceString = strings.at(0);
  int referenceStringLength = referenceString.length();

  for (int j = referenceStringLength; j > 0; --j) {
    string stem = referenceString.substr(0, j);
    bool foundStem = true;
    for (int k = 1; k < numStrings; ++k) {
      if (strings.at(k).find(stem) != 0) {
        foundStem = false;
        break;
      }
    }

    if (foundStem) return stem;
  }

  return "";
}

std::string TabFileValidator::findDecoyPrefix(std::string fileName, int proteinIndex, int labelIndex) {
  std::ifstream file(fileName.c_str());

  if (!file.is_open()) return "error";

  /* Skip header */
  std::string nextRow;
  getline(file, nextRow);

  std::vector<std::string> proteinNames;
  while (getline(file, nextRow)) {
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

  /* Randomly sample 10% of the protein identifiers */
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(proteinNames), std::end(proteinNames), rng);
  int numElements = ceil(0.1 * proteinNames.size());
  proteinNames.resize(numElements);
  
  std::string prefix = getLongestCommonPrefix(proteinNames);
  size_t loc = prefix.find("_");
  if (loc != std::string::npos) {
    prefix = prefix.substr(0, loc + 1);
  }

  if (VERB > 0) {
    std::cerr << "Using protein decoy prefix \"" << prefix << "\"" << std::endl;  
  }
  
  return prefix;
}

std::string TabFileValidator::detectDecoyPrefix(std::string fileName) {
  if (VERB > 0) {
    std::cerr << "Finding protein decoy prefix for " << fileName << std::endl;
  }
  int proteinIndex = -1;
  int labelIndex = -1;

  getProteinAndLabelColumnIndices(fileName, proteinIndex, labelIndex);
  if (proteinIndex == -1 || labelIndex == -1) {
    if (VERB > 0) {
        std::cerr <<  proteinIndex << " " << labelIndex << std::endl;
        std::cerr << "Couldn't find 'Proteins' or 'Label' column in tab-file" << std::endl;
    }
    return "error";
  }
  return findDecoyPrefix(fileName, proteinIndex, labelIndex);
}

bool TabFileValidator::validateTabFiles(std::vector<std::string> files, std::string& decoyPrefix) {
  if (decoyPrefix == "auto") {
    std::string tmpDecoyPrefix = getDecoyPrefix(files);
    if (tmpDecoyPrefix == "error") {
      return false;
    }
    decoyPrefix = tmpDecoyPrefix;
  }
  
  if (!isTabFiles(files)) {
    return false;
  }
  return true;
}