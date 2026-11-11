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

#include <cmath>
#include <cstring>
#include <string_view>

bool TabFileValidator::isTabFile(std::string fileName) {
  std::ifstream file(fileName.c_str());
  
  if (!file.is_open()) {
    if (VERB > 0) {
      std::cerr << "ERROR: Unable to read " << fileName << ", check if the file exists." << std::endl;
    }
    return false;
  }
  
  std::string tmp;
  std::istream &in = std::getline(file, tmp);
  
  if (!in) {
    if (VERB > 0) {
      std::cerr << "ERROR: Unable to read " << fileName << ", check if you have permission to read the file." << std::endl;
    }
    return false;
  }
  
  // try to find a '\t', return true if '\t' is found within the string
  bool isTab = std::find(tmp.begin(), tmp.end(), '\t')!= tmp.end();
  if (!isTab && VERB > 0) {
    std::cerr << "ERROR: " << fileName << " is not tab delimited." << std::endl;
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
  // the prefix is detected in the first file only
  if (fileNames.empty()) {
    return std::string();
  }
  return detectDecoyPrefix(fileNames.front());
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
    for (size_t k = 1; k < numStrings; ++k) {
      if (strings.at(k).find(stem) != 0) {
        foundStem = false;
        break;
      }
    }

    if (foundStem) return stem;
  }

  return "";
}

/*
 * Moves p to the start of the column skip columns further on in a row.
 * Returns false if the row ends before. Like TabReader, a row ends at its
 * first NUL character.
 */
static bool skipColumns(const char*& p, int skip) {
  for (; skip > 0; --skip) {
    p = strchr(p, '\t');
    if (p == NULL) return false;
    ++p;
  }
  return true;
}

// the column that starts at p
static std::string_view column(const char* p) {
  const char* end = strchr(p, '\t');
  return end == NULL ? std::string_view(p)
                     : std::string_view(p, static_cast<size_t>(end - p));
}

std::string TabFileValidator::findDecoyPrefix(std::string fileName, int proteinIndex, int labelIndex) {
  std::ifstream file(fileName.c_str());

  if (!file.is_open()) return "error";

  /* Skip header */
  std::string nextRow;
  getline(file, nextRow);

  // The protein column of the rows with the label -1, one after the other in
  // a single buffer: decoy i is decoyProteins[starts[i], starts[i + 1]).
  // A row only counts as a decoy if its label column comes first.
  std::string decoyProteins;
  std::vector<size_t> starts(1, 0u);
  if (labelIndex < proteinIndex) {
    while (getline(file, nextRow)) {
      const char* p = nextRow.c_str();
      if (!skipColumns(p, labelIndex) || column(p) != "-1" ||
          !skipColumns(p, proteinIndex - labelIndex)) {
        continue;
      }
      decoyProteins.append(column(p));
      starts.push_back(decoyProteins.size());
    }
  }
  const size_t numDecoys = starts.size() - 1;

  /* Randomly sample 10% of the protein identifiers */
  // std::shuffle swaps elements in an order that depends only on the number
  // of elements and the generator, so shuffling the indices selects the same
  // proteins as shuffling the protein names.
  std::vector<size_t> order(numDecoys);
  for (size_t i = 0; i < numDecoys; ++i) order[i] = i;
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(order), std::end(order), rng);
  size_t numElements = static_cast<size_t>(ceil(0.1 * numDecoys));

  // the longest common prefix of the sampled protein names
  std::string prefix;
  if (numElements > 0) {
    const std::string_view reference(decoyProteins.data() + starts[order[0]],
                                     starts[order[0] + 1] - starts[order[0]]);
    size_t length = reference.size();
    for (size_t i = 1; i < numElements && length > 0; ++i) {
      const std::string_view name(decoyProteins.data() + starts[order[i]],
                                  starts[order[i] + 1] - starts[order[i]]);
      size_t common = 0;
      const size_t maxCommon = std::min(length, name.size());
      while (common < maxCommon && name[common] == reference[common]) ++common;
      length = common;
    }
    prefix = std::string(reference.substr(0, length));
  }
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
  if (!isTabFiles(files)) {
    return false;
  }
  
  if (decoyPrefix == "auto") {
    std::string tmpDecoyPrefix = getDecoyPrefix(files);
    if (tmpDecoyPrefix == "error") {
      return false;
    }
    decoyPrefix = tmpDecoyPrefix;
  }
  
  return true;
}