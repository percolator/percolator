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

void TabFileValidator::getProteinAndLabelColumnIndices(std::string file_name, int &proteinIndex, int &labelIndex) {
  // open C++ stream to file
  std::ifstream file(file_name.c_str());

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

std::string TabFileValidator::getLongestCommonPrefix(std::vector<std::string> string_array) {
  size_t num_strings = string_array.size();
  if (num_strings == 0) {
    return "";
  }

  std::string reference_string = string_array.at(0);
  int reference_string_length = reference_string.length();

  for (int j = reference_string_length; j > 0; --j) {
    string stem = reference_string.substr(0, j);
    bool foundStem = true;
    for (int k = 1; k < num_strings; ++k) {
      if (string_array.at(k).find(stem) != 0) {
        foundStem = false;
        break;
      }
    }

    if (foundStem) return stem;
  }

  return "";
}

std::string TabFileValidator::findDecoyPrefix(std::string file_name, int proteinIndex, int labelIndex) {
  std::ifstream file(file_name.c_str());

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

  /* Shuffle protein ids  */
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(proteinNames), std::end(proteinNames), rng);
  /* Check prefix for 10% of all protein ids  */
  int n_elements = ceil(0.1 * proteinNames.size());
  proteinNames.resize(n_elements);
  
  std::string prefix = getLongestCommonPrefix(proteinNames);
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
  int proteinIndex = -1;
  int labelIndex = -1;

  getProteinAndLabelColumnIndices(file_name, proteinIndex, labelIndex);
  if (proteinIndex == -1 || labelIndex == -1) {
    if (VERB > 0) {
        std::cerr <<  proteinIndex << " " << labelIndex << std::endl;
        std::cerr << "Couldn't find 'Proteins' or 'Label' column in tab-file" << std::endl;
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