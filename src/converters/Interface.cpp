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

#include "Interface.h"


Interface::Interface()
{

}

Interface::~Interface()
{

}

string Interface::extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << "on " << host << endl;
  return oss.str();
}

string Interface::greeter() {
  ostringstream oss;
  oss << "\nPin-converter version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2013 Lukas Käll. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukas.kall@scilifelab.se) in the" << endl;
  oss << "School of Biotechnology at KTH - Royal Institute of Technology, Stockholm."
      << endl;
  return oss.str();
}

bool Interface::parseOpt(int argc, char **argv,const std::string &usage) 
{
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl;
  intro << usage;

  // init
  CommandLineParser cmd(intro.str());

  cmd.defineOption("o",
      "outputXML",
      "save output in an XML file",
      "filename");
  cmd.defineOption("m",
      "matches",
      "Maximal number of matches to take in consideration per spectrum",
      "number");
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("a",
      "aa-freq",
      "Calculate amino acid frequency features",
      "",
      TRUE_IF_SET);
  cmd.defineOption("b",
      "PTM",
      "Calculate feature for number of post-translational modifications",
      "",
      TRUE_IF_SET);
  cmd.defineOption("e",
      "enzyme",
      "Type of enzyme \"no_enzyme\",\"elastase\",\"pepsin\",\"proteinasek\",\"thermolysin\",\"chymotrypsin\",\"lys-n\",\"lys-c\",\"arg-c\",\"asp-n\",\"glu-c\",\"trypsin\" default=\"trypsin\"",
      "",
      "trypsin");
  cmd.defineOption("N",
      "PNGaseF",
      "Calculate feature based on N-linked glycosylation pattern resulting from a PNGaseF treatment. (N[*].[ST])",
      "",
      TRUE_IF_SET);
  cmd.defineOption("2",
      "ms2-file",
      "File containing spectra and retention time. The file could be in mzXML, MS2 or compressed MS2 file.",
      "filename");
  cmd.defineOption("M",
      "isotope",
      "Mass difference calculated to closest isotope mass rather than to the average mass.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("p",
      "psm-annotation",
      "An anotation scheme used to convert the psms from the search. An example if Q# was used to describe pyro-glu formation (UNIMOD:28), and S* and T* was used to describe phosphorylation (UNIMOD:21), we would use the option -p *:21:#:28",
      "Scheme");
  cmd.defineOption("P",
      "pattern",
      "Pattern used to identify the decoy PSMs",
      "",
      "pattern");
  /** new parameters for reading the fasta database to obtain the proteins **/
  cmd.defineOption("F",
      "databases",
      "Link to the fasta database/s used in the search against the spectra file/s <target.fasta,[decoy.fasta]> (Including this option will add the proteins to the generated pin file).",
      "",
      "filename");
  cmd.defineOption("c",
      "cleavages",
      "Number of allowed miss cleavages used in the search engine (default 0)(Only valid when using option -F).",
      "",
      "number");
  cmd.defineOption("l",
      "min-length",
      "Minimum peptide length allowed used in the search engine (default 6)(Only valid when using option -F).",
      "",
      "number");
  cmd.defineOption("t",
      "max-length",
      "Maximum peptide length allowed used in the search engine (default 40)(Only valid when using option -F).",
      "",
      "number");
  cmd.defineOption("w",
      "min-mass",
      "Minimum peptide mass allowed used in the search engine (default 400)(Only valid when using option -F).",
      "",
      "number");
  cmd.defineOption("x",
      "max-mass",
      "Maximum peptide mass allowed used in the search engine (default 6000)(Only valid when using option -F).",
      "",
      "number");
  
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  // now query the parsing results

  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (VERB > 0) {
    cerr << extendedGreeter();
  }

  if (cmd.optionSet("o")) {
    xmlOutputFN = cmd.options["o"];
  }
  //option e has been changed, see above
  if (cmd.optionSet("e")) {
    if( cmd.options["e"] == "no_enzyme")
      Enzyme::setEnzyme(Enzyme::NO_ENZYME);
    else if( cmd.options["e"] == "elastase")
      Enzyme::setEnzyme(Enzyme::ELASTASE);
    else if( cmd.options["e"] == "chymotrypsin")
      Enzyme::setEnzyme(Enzyme::CHYMOTRYPSIN);
    else if( cmd.options["e"] == "thermolysin")
      Enzyme::setEnzyme(Enzyme::THERMOLYSIN);
    else if( cmd.options["e"] == "proteinasek")
      Enzyme::setEnzyme(Enzyme::PROTEINASEK);
    else if( cmd.options["e"] == "pepsin")
      Enzyme::setEnzyme(Enzyme::PEPSIN);
    else if( cmd.options["e"] == "trypsin")
      Enzyme::setEnzyme(Enzyme::TRYPSIN);
    else if( cmd.options["e"] == "lys-n")
      Enzyme::setEnzyme(Enzyme::LYSN);
    else if( cmd.options["e"] == "lys-c")
      Enzyme::setEnzyme(Enzyme::LYSC);
    else if( cmd.options["e"] == "arg-c")
      Enzyme::setEnzyme(Enzyme::ARGC);
    else if( cmd.options["e"] == "asp-n")
      Enzyme::setEnzyme(Enzyme::ASPN);
    else if( cmd.options["e"] == "glu-c")
      Enzyme::setEnzyme(Enzyme::GLUC);
    else {
      cerr << "Did not recognize choice of enzyme, uses Trypsin" << endl;
      Enzyme::setEnzyme(Enzyme::TRYPSIN);}
  }
  if (cmd.optionSet("N")) {
    parseOptions.pngasef=true;
  }
  if (cmd.optionSet("a")) {
    parseOptions.calcAAFrequencies=true;
  }
  if (cmd.optionSet("b")) {
    parseOptions.calcPTMs=true;
  }
  if (cmd.optionSet("m")) {
    int m = cmd.getInt("m", 1, 30000);
    parseOptions.hitsPerSpectrum=m;
  }
  if (cmd.optionSet("2")) {
    spectrumFile = cmd.options["2"];
  }
  if (cmd.optionSet("M")) {
    parseOptions.monoisotopic = true;
  }
  if (cmd.optionSet("p")) {
    std::vector<std::string> strs;
    boost::split(strs, cmd.options["p"], boost::is_any_of(":,"));
    if (strs.size()<2) {cerr << "Scheme is malformated" << endl; return 0;}
    for(unsigned int ix=0; ix+1<strs.size(); ix+=2) {
      parseOptions.ptmScheme[strs[ix][0]]=std::stoi(strs[ix+1]);   
      if (VERB > 0) {
        cerr << "Interpreting " << strs[ix][0] << " as modification UNIMOD:" << parseOptions.ptmScheme[strs[ix][0]] << endl; 
      }
    }
  }
  
  if (cmd.optionSet("P")) 
  {
    parseOptions.reversedFeaturePattern = cmd.options["P"];
  }
  
  if (cmd.optionSet("F"))
  {
    //NOTE I do not like this, I should make two parameters, one for target db and another one for decoy db
    std::vector<std::string> strs;
    boost::split(strs, cmd.options["F"], boost::is_any_of(","));
    strs.push_back("");
    parseOptions.targetDb = strs[0];
    parseOptions.decoyDb = strs[1];
    parseOptions.readProteins = true;
  }
  
  if (cmd.optionSet("c"))
  {
    parseOptions.missed_cleavages = cmd.getInt("c", 0, 10);
  }
  
  if (cmd.optionSet("l"))
  {
    parseOptions.peptidelength = cmd.getInt("l",4,20);
  }
  
  if (cmd.optionSet("t"))
  {
    parseOptions.maxpeplength = cmd.getInt("l",6,100);
  }
  
  if (cmd.optionSet("w"))
  {
    parseOptions.minmass = cmd.getInt("l",100,1000);
  }
  
  if (cmd.optionSet("x"))
  {
    parseOptions.maxmass = cmd.getInt("l",100,10000);
  }
  
  if (cmd.arguments.size() > 0)
  {
    targetFN = cmd.arguments[0];
  }
  if (cmd.arguments.size() > 1) 
  {
    decoyFN = cmd.arguments[1];
  }
  if(targetFN == "" && decoyFN == "")
  {
    std::cerr << "Error, one of the input files is missing.\n"; 
    return 0; 
  }
  else if(targetFN != "" && decoyFN == "")
  {
    parseOptions.iscombined = true;
    if(parseOptions.reversedFeaturePattern == "")
    {
      parseOptions.reversedFeaturePattern = "random";
    }
  }
  else
  {
    parseOptions.iscombined = false;
  }
  
  // if there are no arguments left...
  if (cmd.arguments.size() == 0) {
      cerr << "Error: too few arguments.";
      cerr << "\nInvoke with -h option for help\n";
      return 0;
  }
  
  return true;
}

