/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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
#include<string>
#include<map>
#include<vector>
using namespace std;

typedef enum {FALSE_IF_SET=0, TRUE_IF_SET, VALUE, MAYBE} OptionOption;

class Option
{
  public:
    Option(string shrt, string lng, string dest,string hlp="",string hlpType="", OptionOption type=VALUE, string defau="");
    ~Option();
    bool operator == (const string & option);
    OptionOption type;
    string shortOpt;
    string longOpt;
    string help;
    string name;
    string helpType;
    string deflt;
};

class CommandLineParser
{
  public:
    CommandLineParser(string usage="", string tail="");
    ~CommandLineParser();
    void error(string msg);
    void defineOption(string shortOpt, string longOpt, string help="", string helpType="", OptionOption type=VALUE,string defaultVal="");
    void defineOption(string shortOpt, string longOpt, string help, string helpType, string defaultVal) {defineOption(shortOpt,longOpt,help,helpType, VALUE,defaultVal);}
    void parseArgs(int argc, char **argv);
    bool optionSet(string dest) {return (options[dest].length()>0);}
    double getDouble(string dest,double lower,double upper);
    int getInt(string dest,int lower,int upper);
    void help();
    void htmlHelp();
    map<string,string> options;
    vector<string> arguments;
  private:
        unsigned int optMaxLen;
        const static unsigned int lineLen = 80;
    string header, endnote;
    vector<Option> opts;
    void findOption (char **argv, int &index);
};

#endif /*OPTION_H_*/
