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
  
  /* Use these values in place of defineOption()'s shortOpt parameter.
   * NO_SHORT_OPT defines an option with no short flag.
   * EXPERIMENTAL_FEATURE labels the option as experimental in help.
   */
  static const std::string NO_SHORT_OPT, EXPERIMENTAL_FEATURE; 
};

/* The command-line parser.
 */
class CommandLineParser {
 public:
  /* usage and tail define opening and closing strings for help text
   * and errors messages.
   */
  CommandLineParser(std::string usage = "", std::string tail = "");
  ~CommandLineParser();

  /* Fail with an error in the cmdline options.
   */
  void error(std::string msg);

  /* Define a single cmdline option. help is a string describing the
   * option's behavior. helpType, if present, holds a variable name
   * standing for the option's argument. defaultVal describes the
   * default setting if the option is unset.
   */
  void defineOption(std::string shortOpt, std::string longOpt, std::string help = "",
                    std::string helpType = "", OptionOption type = VALUE,
                    std::string defaultVal = "");

  /* If the type parameter is omitted, assume VALUE type. */
  void defineOption(std::string shortOpt, std::string longOpt, std::string help,
                    std::string helpType, std::string defaultVal) {
    defineOption(shortOpt, longOpt, help, helpType, VALUE, defaultVal);
  }

  /* Parse the command-line options. Upon return, the options can be
   * queried via options and/or isOptionSet(). The remaining cmdline
   * arguments will be stored in arguments.
   */
  void parseArgs(int argc, char** argv);

  /* Read the given text file to find more cmdline options.
   */
  void parseArgsParamFile(const std::string paramFile);

  /* Return true if an option is present.
   */
  inline bool isOptionSet(std::string dest) {
    //return (options.find(dest) != options.end());
    return (options[dest].length() > 0);
  }

  /* Return an option's argument as a numerical value. If the argument
   * cannot be so parsed, or if the parsed value is outside of the
   * range given by [lower, upper], an error messages is displayed and
   * the program exits.
   */
  double getDouble(std::string dest, double lower, double upper);
  int getInt(std::string dest, int lower, int upper);
  unsigned int getUInt(std::string dest, int lower, int upper);

  /* Output the help text and exit. */
  void help();

  /* Output the help text as an HTML document and exit. */
  void htmlHelp();

  /* After parsing, use these to look up the options and arguments
   * present on the cmdline. (Boolean options will have "0" or "1" as
   * their value.)
   */
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
