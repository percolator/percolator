#include "msgfdb2Pin.h"

msgfdb2Pin::msgfdb2Pin() {
  
}

msgfdb2Pin::~msgfdb2Pin() {
  if(reader)
    delete reader;
}

string msgfdb2Pin::extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << "on " << host << endl;
  return oss.str();
}

string msgfdb2Pin::greeter() {
  ostringstream oss;
  oss << "\nmsgfdb2pin version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2012 Lukas Käll. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukask@cbr.su.se) in the" << endl;
  oss << "Department of Biochemistry and Biophysics at the Stockholm University."
      << endl;
  return oss.str();
}

bool msgfdb2Pin::parseOpt(int argc, char **argv) {
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   msgfdb2pin [options] -o output.xml target_file decoy_file " << endl << endl;
  intro << "Where output.xml is where the output will be written (ensure to have read and " << endl;
  intro << "write access on the file).target_file is the target msgfdb-file, and decoy_file is" << endl;
  intro << "the decoy msgfdb-file. Small data sets may be merged by replace the msgfdb-files with" << endl;
  intro << "meta files. Meta files are text files containing the paths of msgfdb-files, one" << endl;
  intro << "path per line. For successful result, the different runs should be generated" << endl;
  intro << "under similar condition." << endl;

  // init
  CommandLineParser cmd(intro.str());

  cmd.defineOption("o",
      "outputXML",
      "save output in an XML file",
      "filename");
  //TODO Add an option for which feature should be used for matches per spectra sorting
  cmd.defineOption("m",
      "matches",
      "Maximal number of matches to take in consideration per spectrum, in the order of apperance.",
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
      "Type of enzyme \"no_enzyme\",\"elastase\",\"pepsin\",\"proteinasek\",\"thermolysin\",\"chymotrypsin\",\"trypsin\" default=\"trypsin\"",
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
  /** new parameters for reading the fasta to obtain the proteins **/
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
  if (cmd.optionSet("e")) {
    if( cmd.options["e"] == "no enzyme") 
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
    else  
      Enzyme::setEnzyme(Enzyme::TRYPSIN);
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
    if (strs.size()<2) {cerr << "Scheme is malformated" << endl; exit(-1);}
    for(unsigned int ix=0; ix+1<strs.size(); ix+=2) {
      parseOptions.ptmScheme[strs[ix][0]]=boost::lexical_cast<int>(strs[ix+1]);   
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
    exit(-1); 
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
      exit(-1); // ...error
  }
  
  return true;
}

int msgfdb2Pin::run() {
  // Content of msgfdb files is merged: preparing to write it to xml file
  ofstream xmlOutputStream;
  xmlOutputStream.open(xmlOutputFN.c_str());
  if(!xmlOutputStream && xmlOutputFN != ""){
    cerr << "ERROR: invalid path to output file: " << xmlOutputFN << endl;
    cerr << "Please invoke msgfdb2pin with a valid -o option" << endl;
    exit(-1);
  }
  
  //initialize reader
  parseOptions.targetFN = targetFN;
  parseOptions.decoyFN = decoyFN;
  parseOptions.call = call;
  parseOptions.spectrumFN = spectrumFile;
  parseOptions.xmlOutputFN = xmlOutputFN;
  reader = new msgfdbReader(&parseOptions);
  
  reader->init();
  reader->print(xmlOutputStream);
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return 0;
}

int main(int argc, char** argv) {
  msgfdb2Pin* pmsgfdb2Pin = new msgfdb2Pin();
  int retVal = -1;
  if (pmsgfdb2Pin->parseOpt(argc, argv)) {
    retVal = pmsgfdb2Pin->run();
  }
  delete pmsgfdb2Pin;
  Globals::clean();
  return retVal;
}