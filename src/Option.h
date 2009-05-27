/*******************************************************************************
 Copyright (c) 2008-9 Lukas Käll

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: Option.h,v 1.11 2009/05/27 07:24:08 lukall Exp $
 
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
