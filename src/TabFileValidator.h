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
    static bool isTabFile(std::string fileName);
    static bool isTabFiles(std::vector<std::string> fileNames);
    static std::string getDecoyPrefix(std::vector<std::string> fileNames);
    static std::string detectDecoyPrefix(std::string fileName);
    static bool validateTabFiles(std::vector<std::string> fileNames, std::string &decoyPrefix);
    static void getProteinAndLabelColumnIndices(std::string fileName, int &proteinIndex, int &labelIndex);
    static std::string findDecoyPrefix(std::string fileName, int proteinIndex, int labelIndex);
    static std::string getLongestCommonPrefix(std::vector<std::string> strings);
};

#endif /*TABFILEVALIDATOR_H_*/