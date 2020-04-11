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
#ifndef OPTION_H_
#define OPTION_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <cctype>

#include "MyException.h"

typedef enum {
  FALSE_IF_SET = 0, TRUE_IF_SET, VALUE, MAYBE
} OptionOption;

class Option {
 public:
  Option(std::string shrt, std::string lng, std::string dest, 
         std::string hlp = "", std::string hlpType = "", 
         OptionOption type = VALUE, std::string defau = "");
  ~Option();
  bool operator ==(const std::string& option);
  OptionOption type;
  std::string shortOpt;
  std::string longOpt;
  std::string help;
  std::string name;
  std::string helpType;
  std::string deflt;
  
  /* allow options without a short flag and for experimental features */
  static const std::string NO_SHORT_OPT, EXPERIMENTAL_FEATURE; 
};

class CommandLineParser {
 public:
  CommandLineParser(std::string usage = "", std::string tail = "");
  ~CommandLineParser();
  void error(std::string msg);
  void defineOption(std::string shortOpt, std::string longOpt, std::string help = "",
                    std::string helpType = "", OptionOption type = VALUE,
                    std::string defaultVal = "");
  void defineOption(std::string shortOpt, std::string longOpt, std::string help,
                    std::string helpType, std::string defaultVal) {
    defineOption(shortOpt, longOpt, help, helpType, VALUE, defaultVal);
  }
  void parseArgs(int argc, char** argv);
  void parseArgsParamFile(const std::string paramFile);
  inline bool optionSet(std::string dest) {
    //return (options.find(dest) != options.end());
    return (options[dest].length() > 0);
  }
  double getDouble(std::string dest, double lower, double upper);
  int getInt(std::string dest, int lower, int upper);
  int getIntLowerBound(string dest, int lower);
  void help();
  void htmlHelp();
  std::map<std::string, std::string> options;
  std::vector<std::string> arguments;
 private:
  size_t optMaxLen;
  const static unsigned int lineLen = 80;
  std::string header, endnote;
  std::vector<Option> opts;
  void findOption(char** argv, int& index, int argc);
  static inline std::string &rtrim(std::string &s);
};

#endif /*OPTION_H_*/
