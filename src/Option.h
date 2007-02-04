/*******************************************************************************
 * Percolator v 1.00
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas KŠll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Option.h,v 1.2 2007/02/04 04:33:21 lukall Exp $
 *******************************************************************************/
#ifndef OPTION_H_
#define OPTION_H_

typedef enum {FALSE_IF_SET=0, TRUE_IF_SET, VALUE} OptionOption;

class Option
{
  public:
    Option(string shrt, string lng, string dest,string hlp="",string hlpType="", OptionOption type=VALUE);
    ~Option();
    bool operator == (const string & option);
    OptionOption type;
    string shortOpt;
    string longOpt;
    string help;    
    string name;
    string helpType;
};

class CommandLineParser
{
  public:
    CommandLineParser(string usage="");
    ~CommandLineParser();
    void error(string msg);
    void defineOption(string shortOpt, string longOpt, string help="", string helpType="", OptionOption type=VALUE,string defaultVal="");
    void defineOption(string shortOpt, string longOpt, string help, string helpType, string defaultVal) {defineOption(shortOpt,longOpt,help,helpType, VALUE,defaultVal);}
    void parseArgs(int argc, char **argv);
    bool optionSet(string dest) {return (options[dest].length()>0);}
    double getDouble(string dest,double lower,double upper);
    int getInt(string dest,int lower,int upper);
    void help();
    map<string,string> options;
    vector<string> arguments;
  private:
        unsigned int optMaxLen; 
        const static unsigned int lineLen = 80;
    string header;      
    vector<Option> opts;
    void findOption (char **argv, int &index);
};

#endif /*OPTION_H_*/
