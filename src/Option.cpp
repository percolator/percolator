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
  if (strlen(longOpt.c_str())+strlen(shortOpt.c_str())+strlen(helpType.c_str())> optMaxLen)
    optMaxLen = strlen(longOpt.c_str())+strlen(shortOpt.c_str())+strlen(helpType.c_str());
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
  cerr << header << endl << "Options:" << endl;
  int descLen = optMaxLen + 12;
  int helpLen = lineLen - descLen;
    for (unsigned int i=0; i < opts.size(); i++) {
      unsigned int j=0;
      string desc = " " + opts[i].shortOpt + " or " + opts[i].longOpt;
      if (opts[i].helpType.length()>0) {
        desc += " <" + opts[i].helpType + ">";
      }
      cerr.width(descLen);
      cerr << left << desc; 
      cerr.width(0);
      cerr << opts[i].help.substr(j,helpLen) << endl; 
      while ((j+=helpLen)<opts[i].help.length()) {
        cerr.width(descLen);
        cerr << " "; 
        cerr.width(0);
        cerr << opts[i].help.substr(j,helpLen) << endl; 
      }
  }
  exit(0);
}

void CommandLineParser::findOption(char **argv, int &index) {
  if ((string)argv[index] == "-h" || (string)argv[index] == "--help")
    help();

  string optstr = (string)argv[index];
  string valstr("");
  unsigned int eqsign=optstr.find('=');
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
