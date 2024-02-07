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
#ifndef TABFILEVALIDATOR_H_
#define TABFILEVALIDATOR_H_

#include <string>
#include <vector>
#include <random>
#include <algorithm>

#include "Globals.h"
#include "TabReader.h"

class TabFileValidator {
  public:
    static bool isTabFile(std::string file_name);
    static bool isTabFiles(std::vector<std::string> files);
    std::string getDecoyPrefix(std::vector<std::string> fileList);
    std::string detectDecoyPrefix(std::string file_name);
    bool validateTabFiles(std::vector<std::string> files, std::string* decoy_prefix);
    void getProteinAndLabelColumnIndices(std::string file_name, int &proteinIndex,int &labelIndex);
    std::string findDecoyPrefix(std::string file_name, int proteinIndex, int labelIndex);
    std::string getLongestCommonPrefix(std::vector<std::string> arr);
    static bool decoyWarningTripped;
};

#endif /*TABFILEVALIDATOR_H_*/