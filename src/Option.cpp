/*******************************************************************************
 * Percolator v 1.00
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Option.cpp,v 1.5 2007/02/04 04:50:39 lukall Exp $
 *******************************************************************************/
#include <iostream>
#include <string>
#include <map>
#include <vector>
using namespace std;
#include "Option.h"

Option::Option (string shrt, string lng, string nm, string hlp,string hlpType, OptionOption typ) {
  type=typ;
  shortOpt=shrt;
  longOpt=lng;
  help=hlp;
  helpType=hlpType;
  name=nm;
}

Option::~Option () {}

bool Option::operator == (const string & option) {
  return (shortOpt == option || longOpt == option);
}

CommandLineParser::CommandLineParser(string usage) {
  header=usage;
  optMaxLen=0;
  defineOption("h","help","Display this message");
}

CommandLineParser::~CommandLineParser() {}

double CommandLineParser::getDouble(string dest,double lower,double upper){
    double val = atof(options[dest].c_str());
    if (val < lower || val > upper) {
      cerr << "-" << dest << " option requres a float between " << lower << " and " << upper << endl;
      exit(-1); 
    }
    return val;
}

int CommandLineParser::getInt(string dest,int lower,int upper) {
    int val = atoi(options[dest].c_str());
    if (val < lower || val > upper) {
      cerr << "-" << dest << " option requres an integer between " << lower << " and " << upper << endl;
      exit(-1); 
    }
    return val;
}


void CommandLineParser::defineOption (string shortOpt, string longOpt, string help, string helpType, OptionOption typ, string dfault) {
  opts.insert(opts.begin(),Option("-"+shortOpt, "--" + longOpt, shortOpt, help,helpType, typ));
  options[shortOpt] = dfault;
  if (longOpt.length()+helpType.length()> optMaxLen)
    optMaxLen=longOpt.length()+helpType.length();
}

void CommandLineParser::parseArgs(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-')
      findOption(argv, i);
    else
      arguments.insert(arguments.end(), argv[i]);
  }
}

void CommandLineParser::error(string msg) {
  cerr << header << endl << msg << endl;
  exit(1);
}

void CommandLineParser::help() {
  string::size_type descLen = optMaxLen + 8;
  string::size_type helpLen = lineLen - descLen;
  cerr << header << endl << "Options:" << endl;
  for (unsigned int i=opts.size(); i-- ;) {
    string::size_type j=0;
    cerr << " " << opts[i].shortOpt;
    if (opts[i].helpType.length()>0)
      cerr << " <" << opts[i].helpType << ">";
    cerr << endl;
    string desc = " " + opts[i].longOpt;
    if (opts[i].helpType.length()>0) {
      desc += " <" + opts[i].helpType + ">";
    }
    while (j<opts[i].help.length()) {
      cerr.width(descLen);
      cerr << left << desc;
      desc = " "; 
      cerr.width(0);
      string::size_type l = helpLen;
      if (j+l<opts[i].help.length()) {
        string::size_type p = opts[i].help.rfind(' ',j+l);
        if (p != string::npos && p>j) 
          l = p-j+1;
      }
      cerr << opts[i].help.substr(j,l) << endl;
      j += l; 
    }
  }
  exit(0);
}

void CommandLineParser::findOption(char **argv, int &index) {
  if ((string)argv[index] == "-h" || (string)argv[index] == "--help")
    help();

  string optstr = (string)argv[index];
  string valstr("");
  string::size_type eqsign=optstr.find('=');
  if (eqsign!=string::npos) {
    valstr = optstr.substr(eqsign+1);
    optstr = optstr.substr(0,eqsign);
  }    
  for (unsigned int i = 0; i < opts.size(); i++) {
    if (opts[i] == optstr) {
      switch (opts[i].type) {
        case FALSE_IF_SET:
          options[opts[i].name] = "0";
          break;
        case TRUE_IF_SET:
          options[opts[i].name] = "1";
          break;
        case VALUE:
          if (valstr.length()>0) {
            options[opts[i].name] = valstr;        
          } else {
            options[opts[i].name] = argv[index+1];
            index++;
          }
          break;
        default:
          break;
      };
      return;
    }
  }
  error("Error: invalid argument");
}
