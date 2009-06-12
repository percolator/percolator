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

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
#include "Option.h"

template <class T>
bool from_string(T& t,
                 const std::string& s)
{
  std::istringstream iss(s);
  return !(iss >> t).fail();
}


void searchandreplace( string& source, const string& find, const string& replace )
{
  size_t j;
  for (;(j = source.find( find )) != string::npos;)
    source.replace( j, find.length(), replace );
}

Option::Option (string shrt, string lng, string nm, string hlp,string hlpType, OptionOption typ, string dfl) {
  type=typ;
  shortOpt=shrt;
  longOpt=lng;
  help=hlp;
  helpType=hlpType;
  name=nm;
  deflt = dfl;
}

Option::~Option () {}

bool Option::operator == (const string & option) {
  return (shortOpt == option || longOpt == option);
}

CommandLineParser::CommandLineParser(string usage, string tail) {
  header=usage;
  endnote=tail;
  optMaxLen=0;
  defineOption("h","help","Display this message");
}

CommandLineParser::~CommandLineParser() {}

double CommandLineParser::getDouble(string dest,double lower,double upper){
    double val;
    from_string<double>(val, options[dest]);
    if (!from_string<double>(val, options[dest]) || (val < lower || val > upper)) {
      cerr << "-" << dest << " option requires a float between " << lower << " and " << upper << endl;
      exit(-1);
    }
    return val;
}

int CommandLineParser::getInt(string dest,int lower,int upper) {
    int val;
    if (!from_string<int>(val, options[dest]) || val < lower || val > upper) {
      cerr << "-" << dest << " option requires an integer between " << lower << " and " << upper << endl;
      exit(-1);
    }
    return val;
}


void CommandLineParser::defineOption (string shortOpt, string longOpt, string help, string helpType, OptionOption typ, string dfault) {
  opts.insert(opts.begin(),Option("-"+shortOpt, "--" + longOpt, shortOpt, help,helpType, typ,dfault));
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
  cerr << endl << endnote << endl;
  exit(0);
}


void CommandLineParser::htmlHelp() {
  cerr << "<html><title>Title</title><body><blockquote>" << endl;
  cerr << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" />" << endl;

  string htmlHeader = header;
  searchandreplace(htmlHeader,"\n","<br/>");
  cerr << htmlHeader << endl << "Options:" << endl;
  cerr << "<table border=0>" << endl;
  for (unsigned int i=opts.size(); i-- ;) {
    cerr << "<tr><td><code>" << opts[i].shortOpt;
    if (opts[i].helpType.length()>0)
      cerr << " &lt;" << opts[i].helpType << "&gt;";
    cerr << "</code>, <code>";
    cerr << " " + opts[i].longOpt;
    if (opts[i].helpType.length()>0) {
      cerr << " &lt;" << opts[i].helpType << "&gt;";
    }
    cerr << "</code></td>" << endl;
    cerr << "<td>" << opts[i].help << "</td></tr>" << endl;
  }
  cerr << "</table>" << endl;
  string htmlEnd = endnote;
  searchandreplace(htmlEnd,"\n","<br>");
  cerr << "<br/>" << endl << htmlEnd << "<br/>" << endl;
  cerr << "</blockquote></body></html>" << endl;
  exit(0);
}

void CommandLineParser::findOption(char **argv, int &index) {
  if ((string)argv[index] == "-html" || (string)argv[index] == "--html")
    htmlHelp();

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
        case MAYBE:
          if (valstr.length()>0) {
            options[opts[i].name] = valstr;
          } else if (argv[index+1][0]!='-') {
            options[opts[i].name] = argv[index+1];
            index++;
          } else {
            options[opts[i].name] = opts[i].deflt;
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
