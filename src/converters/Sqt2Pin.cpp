/*
 * Sqt2Pin.cpp
 *
 *  Created on: Jun 16, 2010
 *      Author: lukask
 */

using namespace std;
#include <iostream>
#include <fstream>
#include "Enzyme.h"
#include "Sqt2Pin.h"
#include "config.h"
#include "serializer.hxx"
#include "MSReader.h"
#include "Spectrum.h"
#include "MSToolkitTypes.h"
#include "MassHandler.h"
#include "SqtReader.h"
#include "DataSet.h"

Sqt2Pin::Sqt2Pin() {
  // TODO Auto-generated constructor stub
}

string Sqt2Pin::greeter() {
  ostringstream oss;
  oss << "Sqt2Pin version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss
  << "Copyright (c) 2010 Lukas Käll. All rights reserved."
  << endl;
  oss << "Written by Lukas Käll (lukask@cbr.su.se) in the" << endl;
  oss << "Department of Biochemistry and Biophysics at the Stockholm University."
      << endl;
  return oss.str();
}


bool Sqt2Pin::parseOpt(int argc, char **argv) {
  xmlOutputFN="";
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   sqt2pin [options] target.sqt decoy.sqt" << endl;
  intro << "   where target.sqt is the target sqt-file," << endl;
  intro << "     and decoy.sqt is the decoy sqt-file," << endl;
  intro
  << "Small data sets may be merged by replace the sqt-files with meta"
  << endl;
  intro
  << "files. Meta files are text files containing the paths of sqt-files, one path"
  << endl;
  intro
  << "per line. For successful result, the different runs should be generated under"
  << endl;
  intro << "similair condition." << endl;
  // init
  CommandLineParser cmd(intro.str());

  cmd.defineOption("o",
      "outputXML",
      "save output in an XML file",
      "filename");
  cmd.defineOption("m",
      "matches",
      "Maximal number of matches to take in consideration per spectrum when using sqt-files",
      "number");
  cmd.defineOption("g",
      "gist-in",
      "Input files are given as gist files. In this case first argument should be a file name \
      of the data file, the second the label file. Labels are interpreted as 1 -- positive train \
      and test set, -1 -- negative train set, -2 -- negative in test set.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("u",
      "unitnorm",
      "Use unit normalization [0-1] instead of standard deviation normalization",
      "",
      TRUE_IF_SET);
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
  cmd.defineOption("Q",
      "quadratic",
      "Calculate quadratic feature terms",
      "",
      TRUE_IF_SET);
  cmd.defineOption("y",
      "notryptic",
      "Turn off calculation of tryptic/chymo-tryptic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("c",
      "chymo",
      "Replace tryptic features with chymo-tryptic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("e",
      "elastase",
      "Replace tryptic features with elastase features.",
      "",
      TRUE_IF_SET);
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

  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  // now query the parsing results

  if (cmd.optionSet("Y")) { tokyoCabinetTmpFN = cmd.options["Y"];
  } else  { tokyoCabinetTmpFN = "/tmp/percolator-tmp.tcb"; }

  if (cmd.optionSet("o")) {
    xmlOutputFN = cmd.options["o"];
  }
  if (cmd.optionSet("Q")) {
    parseOptions.calcQuadraticFeatures=true;
  }
  if (cmd.optionSet("y")) {
    Enzyme::setEnzyme(Enzyme::NO_ENZYME);
  }
  if (cmd.optionSet("e")) {
    Enzyme::setEnzyme(Enzyme::ELASTASE);
  }
  if (cmd.optionSet("c")) {
    Enzyme::setEnzyme(Enzyme::CHYMOTRYPSIN);
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
    MassHandler::setMonoisotopicMass(true);
  }
  if (cmd.arguments.size() > 0) targetFN = cmd.arguments[0];
  if (cmd.arguments.size() > 1) decoyFN = cmd.arguments[1];
  return true;
}

int Sqt2Pin::run() {
  //if ( xmlOutputFN.size() != 0 ) {
    // content of sqt files is merged. Results are printed on both cout and file
    ofstream xmlOutputStream; // stream for output XML file
    xmlOutputStream.open(xmlOutputFN.c_str());

    xercesc::XMLPlatformUtils::Initialize ();
    std::auto_ptr<percolatorInNs::featureDescriptions> fdes_p ( new ::percolatorInNs::featureDescriptions());
    std::auto_ptr< ::percolatorInNs::experiment > ex_p ( new ::percolatorInNs::experiment( "mitt enzym" , fdes_p ));

    int maxCharge = -1;
    int minCharge = 10000;
    FragSpectrumScanDatabase database;

    /*
	   The function "tcbdbopen" in Tokyo Cabinet does not have O_EXCL as is possible in the unix system call open (see "man 2 open").
	   This may be a security issue if the filename to the tokyo cabinet database is in a directory
	   that other users have write access to. They could add a symbolic link pointing somewhere else.

	   It would be better if Tokyo Cabinet would fail if the database existed in our case when we use a tempory file.
     */

    database.init(tokyoCabinetTmpFN);

    if (targetFN != "" && parseOptions.reversedFeaturePattern.empty() ) {
      // First we only search for the maxCharge and minCharge. This done by passing the argument justSearchMaxMinCharge
      SqtReader::translateSqtFileToXML( targetFN,ex_p->featureDescriptions(),  ex_p->fragSpectrumScan(), false /* is_decoy */, parseOptions, &maxCharge, &minCharge, SqtReader::justSearchMaxMinCharge ,  database );
      SqtReader::translateSqtFileToXML( decoyFN, ex_p->featureDescriptions(),  ex_p->fragSpectrumScan(),  true /* is_decoy */, parseOptions, &maxCharge, &minCharge,  SqtReader::justSearchMaxMinCharge , database );
      // Now we do full parsing of the Sqt file, and translating it to XML
      SqtReader::translateSqtFileToXML( targetFN,ex_p->featureDescriptions(),  ex_p->fragSpectrumScan(),  false /* is_decoy */ , parseOptions, &maxCharge, &minCharge,  SqtReader::fullParsing, database  );
      SqtReader::translateSqtFileToXML( decoyFN, ex_p->featureDescriptions(),  ex_p->fragSpectrumScan(),  true /* is_decoy */, parseOptions, &maxCharge, &minCharge,  SqtReader::fullParsing, database  );

    } else {
      // First we only search for the maxCharge and minCharge. This done by passing the argument justSearchMaxMinCharge
      SqtReader::translateSqtFileToXML( targetFN,ex_p->featureDescriptions(),     ex_p->fragSpectrumScan() ,  false /* is_decoy */, parseOptions, &maxCharge, &minCharge, SqtReader::justSearchMaxMinCharge, database );
      // Now we do full parsing of the Sqt file, and translating it to XML
      SqtReader::translateSqtFileToXML( targetFN,ex_p->featureDescriptions(),     ex_p->fragSpectrumScan() ,  true /* is_decoy */, parseOptions, &maxCharge, &minCharge, SqtReader::fullParsing, database );
    }
    //    pCheck = new SqtSanityCheck();
    //    assert(pCheck);

    std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    xmlOutputStream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    std::cout << "<experiment  xmlns=\"" << PERCOLATOR_IN_NAMESPACE <<
        "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\""
        << PERCOLATOR_IN_NAMESPACE <<
        " http://github.com/percolator/percolator/raw/master/src/xml/percolator_in.xsd\">"
        << std::endl;
    xmlOutputStream << "<experiment  xmlns=\"" << PERCOLATOR_IN_NAMESPACE <<
        "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\""
        << PERCOLATOR_IN_NAMESPACE <<
        " http://github.com/percolator/percolator/raw/master/src/xml/percolator_in.xsd\">"
        << std::endl;

    std::string enzyme;
    switch ( Enzyme::getEnzymeType() ) {
    case Enzyme::NO_ENZYME : { enzyme = "NO_ENZYME"; break; }
    case Enzyme::TRYPSIN : { enzyme = "trypsin"; break; }
    case Enzyme::CHYMOTRYPSIN : { enzyme = "chymotrypsin"; break; }
    case Enzyme::ELASTASE : { enzyme = "elastase"; break; }
    }
    std::cout << "   <enzyme>" << enzyme << "</enzyme>" << std::endl;
    xmlOutputStream << "   <enzyme>" << enzyme << "</enzyme>" << std::endl;
    serializer ser;
    ser.start (std::cout);
    ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",
        ex_p->featureDescriptions() );
    database.print(ser);
    std::cout << "</experiment>" << std::endl;
    serializer serXML;
    serXML.start (xmlOutputStream);
    serXML.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",
        ex_p->featureDescriptions() );
    database.print(serXML);
    xmlOutputStream << "</experiment>" << std::endl;
    xmlOutputStream.close(); // close stream for output XML file
  //}
  if (spectrumFile.size() > 0) {
    readRetentionTime(spectrumFile);
  }
  return 0;
}

void Sqt2Pin::readRetentionTime(string filename) {
  MSReader r;
  Spectrum s;
  r.setFilter(MS2);
  char* cstr = new char[filename.size() + 1];
  strcpy(cstr, filename.c_str());
  r.readFile(cstr, s);
  while (s.getScanNumber() != 0) {
    scan2rt[s.getScanNumber()] = (double)s.getRTime();
    r.readFile(NULL, s);
  }
  delete[] cstr;
}


Sqt2Pin::~Sqt2Pin() {
  // TODO Auto-generated destructor stub
}

int main(int argc, char** argv) {
  Sqt2Pin* pSqt2Pin = new Sqt2Pin();
  int retVal = -1;
  if (pSqt2Pin->parseOpt(argc, argv)) {
    retVal = pSqt2Pin->run();
  }
  delete pSqt2Pin;
  Globals::clean();
  return retVal;
}

