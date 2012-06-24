#include "mzidentml2pin.h"
#include <boost/lexical_cast.hpp>


enum enzyme_type { enzyme_type_arg_no_enzyme, enzyme_type_arg_elastase, 
		    enzyme_type_arg_pepsin, enzyme_type_arg_proteinasek, 
		    enzyme_type_arg_thermolysin, enzyme_type_arg_chymotrypsin, 
		    enzyme_type_arg_trypsin};
		    

Mzidentml2pin::Mzidentml2pin()
{

}

Mzidentml2pin::~Mzidentml2pin()
{

}


std::string Mzidentml2pin::greeter()
{
  ostringstream oss;
  oss << "mzidentml2pin version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2012 Lukas Käll. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukask@cbr.su.se) in the" << endl;
  oss << "Department of Biochemistry and Biophysics at the Stockholm University."
      << endl;
  return oss.str();
  
}

std::string Mzidentml2pin::extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << "on " << host << endl;
  return oss.str();
}




bool Mzidentml2pin::parseOpt(int argc, char **argv)
{
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   mzidentml2pin [options] target.sqt decoy.sqt" << endl << endl;
  intro << "Target.sqt is the target sqt-file, and decoy.sqt is" << endl;
  intro << "the decoy sqt-file. Small data sets may be merged by replace the sqt-files with" << endl;
  intro << "meta files. Meta files are text files containing the paths of sqt-files, one" << endl;
  intro << "path per line. For successful result, the different runs should be generated" << endl;
  intro << "under similar condition." << endl;

  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("e",
      "enzyme",
      "Type of enzyme \"no_enzyme\",\"elastase\",\"pepsin\",\"proteinasek\",\"thermolysin\",\"chymotrypsin\",\"trypsin\" default=\"trypsin\"",
      "",
      "trypsin");
  cmd.defineOption("b",
      "ptm",
      "Calculate feature for number of post-translational modifications",
      "",
      TRUE_IF_SET);
  cmd.defineOption("a",
      "aa-freq",
      "Calculate amino acid frequency features",
      "",
      TRUE_IF_SET);
  cmd.defineOption("N",
      "pngasef",
      "Calculate feature based on N-linked glycosylation pattern resulting from a PNGaseF treatment. (N[*].[ST])",
      "",
      TRUE_IF_SET);
  cmd.defineOption("Z",
      "combined",
      "The input will be a target/decoy combined file.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("P",
      "pattern",
      "Pattern used to identify the decoy PSMs",
      "",
      "pattern");
  cmd.defineOption("o",
      "outputXML",
      "save output in an XML file",
      "filename");
  cmd.defineOption("m",
      "matches",
      "Maximal number of matches to take in consideration per spectrum when using sqt-files",
      "number");
  cmd.defineOption("M",
      "isotope",
      "Mass difference calculated to closest isotope mass rather than to the average mass.",
      "",
      TRUE_IF_SET);
  
  cmd.parseArgs(argc, argv);
  
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (VERB > 0) {
    cerr << extendedGreeter();
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
  
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (VERB > 0) {
    cerr << extendedGreeter();
  }

  if (cmd.optionSet("o")) {
    xmlOutputFN = cmd.options["o"];
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

  if (cmd.optionSet("M")) {
    MassHandler::setMonoisotopicMass(true);
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
  
  std::string pattern = "";
  if (cmd.optionSet("P")) pattern = cmd.options["P"];
  bool iscombined = cmd.optionSet("Z");
  
  if(iscombined)
  {
    if(cmd.arguments.size() > 1)
    {
      std::cerr << "Error, there should be only one input file.\n"; 
      exit(-1);
    }
    else if (pattern != "")
    {
      parseOptions.reversedFeaturePattern = pattern;
      parseOptions.iscombined = true;
      targetFN = cmd.arguments[0];
    }
    else
    {
      std::cerr << "Error, pattern should contain a valid set of alphanumberic characters.\n"; 
      exit(-1);
    }
  }
  else
  {
    if (cmd.arguments.size() > 0) targetFN = cmd.arguments[0];
    if (cmd.arguments.size() > 1) decoyFN = cmd.arguments[1];
    
    if(targetFN == "" || decoyFN == "")
    {
      std::cerr << "Error, one of the input files is missing.\n"; 
      exit(-1); 
    }
  }
  // if there are no arguments left...
  if (cmd.arguments.size() == 0) {
      cerr << "Error: too few arguments.";
      cerr << "\nInvoke with -h option for help\n";
      exit(-1); // ...error
  }
  
  
  return true;
}


int Mzidentml2pin::run() {
  // Content of sqt files is merged: preparing to write it to xml file
  ofstream xmlOutputStream;
  xmlOutputStream.open(xmlOutputFN.c_str());
  if(!xmlOutputStream && xmlOutputFN != ""){
    cerr << "ERROR: invalid path to output file: " << xmlOutputFN << endl;
    cerr << "Please invoke mzidentml2pin with a valid -o option" << endl;
    exit(-1);
  }
  
    //initialize reader
  parseOptions.boost_serialization = boost_serialization;
  parseOptions.peptidelength = 6;
  reader = new MzidentmlReader(parseOptions);
  
  xercesc::XMLPlatformUtils::Initialize ();
  // initializing xercesc objects corresponding to pin element...
  // ... <featureDescriptions>
  std::auto_ptr<percolatorInNs::featureDescriptions>
  fdes_p (new ::percolatorInNs::featureDescriptions());

  // ... <process_info>
  percolatorInNs::process_info::command_line_type command_line = call;
  std::auto_ptr<percolatorInNs::process_info>
  proc_info (new ::percolatorInNs::process_info(command_line));

  // ... <experiment>
  std::auto_ptr< ::percolatorInNs::experiment >
  ex_p (new ::percolatorInNs::experiment("mitt enzym", proc_info, fdes_p));

  //NOTE here is where I declare the serializer depending on the library we built converters with
  vector<serialize_scheme*> databases;

  if (!parseOptions.iscombined) 
  {
    
    std::cerr << "Reading input from mzid files:\n";
    
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */,databases,0);
    reader->translateFileToXML(decoyFN, ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), true /* is_decoy */,databases,0);
    
  } 
  else 
  {
    
    std::cerr << "Reading input from a combined (target-decoy) mzid file .." << std::endl;
 
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */,databases,0);
    
  }


  xercesc::XMLPlatformUtils::Terminate();

  string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
  string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
  string headerStr = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n" +
      string("<experiment xmlns=\"") + PERCOLATOR_IN_NAMESPACE + "\"" +
      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" +
      " xsi:schemaLocation=\"" + PERCOLATOR_IN_NAMESPACE +
      " https://github.com/percolator/percolator/raw/pin-" + schema_major +
      "-" + schema_minor + "/src/xml/percolator_in.xsd\"> \n";
  if (xmlOutputFN == "") cout << headerStr;
  else {
    xmlOutputStream << headerStr;
    cerr <<  "The output will be written to " << xmlOutputFN << endl;
  }

  string enzymeStr = "\n<enzyme>" + Enzyme::getStringEnzyme() + "</enzyme>\n";
  if (xmlOutputFN == "") cout << enzymeStr;
  else xmlOutputStream << enzymeStr;

  string commandLine = "\n<process_info>\n" +
      string("  <command_line>") + call.substr(0,call.length()-1)
      + "</command_line>\n" + "</process_info>\n";
  if (xmlOutputFN == "") cout << commandLine;
  else xmlOutputStream << commandLine;

  xercesc::XMLPlatformUtils::Initialize ();

  cerr << "\nWriting output:\n";
  // print to cout (or populate xml file)
  // print features
  {
    serializer ser;
    if (xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",
        ex_p->featureDescriptions());
  }

  // print fragSpecturmScans
  std::cerr << "Databases : " << databases.size() << std::endl;
  for(int i=0; i<databases.size();i++) {
    serializer ser;
    if (xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    if(VERB>1){
      cerr << "outputting content of " << databases[i]->id
          << " (and correspondent decoy file)\n";
    }
    databases[i]->print(ser);
    databases[i]->terminte();
  }

  // print closing tag
  if (xmlOutputFN == "") std::cout << "</experiment>" << std::endl;
  else {
    xmlOutputStream << "</experiment>" << std::endl;
    xmlOutputStream.close();
  }

  xercesc::XMLPlatformUtils::Terminate();

  cerr << "\nAll the input files have been successfully processed"<< endl;

  return 0;
}

int main(int argc, char** argv) {
  Mzidentml2pin* pSqt2Pin = new Mzidentml2pin();
  int retVal = -1;
  if (pSqt2Pin->parseOpt(argc, argv)) {
    retVal = pSqt2Pin->run();
  }
  delete pSqt2Pin;
  Globals::clean();
  return retVal;
}