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

#include "XMLInterface.h"

XMLInterface::XMLInterface() : xmlInputFN(""), schemaValidation(false), otherCall(""), xmlOutputFN(""), reportUniquePeptides(false) {}

int XMLInterface::readPin(SetHandler & setHandler, SanityCheck *& pCheck, ProteinProbEstimator * protEstimator) {    
#ifdef XML_SUPPORT  
  DataSet * targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(1);
  DataSet * decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(-1);

  xercesc::XMLPlatformUtils::Initialize();
  
  try {
    
    using namespace xercesc;
    
    std::ifstream xmlInStream;
    xmlInStream.exceptions(ifstream::badbit | ifstream::failbit);
    xmlInStream.open(xmlInputFN.c_str());

    string schemaDefinition= Globals::getInstance()->getXMLDir()+PIN_SCHEMA_LOCATION+string("percolator_in.xsd");
    parser p;
    xml_schema::dom::auto_ptr<DOMDocument> doc(p.start(
        xmlInStream, xmlInputFN.c_str(), schemaValidation,
        schemaDefinition, PIN_VERSION_MAJOR, PIN_VERSION_MINOR));

    doc = p.next();
    // read enzyme element
    // the enzyme element is a subelement but CodeSynthesis Xsd does not
    // generate a class for it. (I am trying to find a command line option
    // that overrides this decision). As for now special treatment is needed
    char* value = XMLString::transcode(doc->getDocumentElement()->getTextContent());
    
    if(VERB > 1) std::cerr << "enzyme=" << value << std::endl;
    
    Enzyme::setEnzyme(value);
    XMLString::release(&value);
    doc = p.next();

    //checking if database is present to jump it
    bool hasProteins = false;
    if(XMLString::equals(databasesStr, doc->getDocumentElement()->getTagName())) {
      //NOTE I dont really need this info, do I? good to have it though
      /*
        std::unique_ptr< ::percolatorInNs::databases > 
      databases( new ::percolatorInNs::databases(*doc->getDocumentElement()));
      */
      doc = p.next();
      hasProteins = true;
    }
    
    // read process_info element
    percolatorInNs::process_info
    processInfo(*doc->getDocumentElement());
    otherCall = processInfo.command_line();
    doc = p.next();


    if (XMLString::equals(calibrationStr,doc->getDocumentElement()->getTagName())) {
      //NOTE the calibration should define the initial direction
      //percolatorInNs::calibration calibration(*doc->getDocumentElement());
      doc = p.next();
    };

    percolatorInNs::featureDescriptions featureDescriptions(*doc->getDocumentElement());

    //I want to get the initial values that are present in feature descriptions
    std::vector<double> init_values;
    for (const auto & descr : featureDescriptions.featureDescription()) {
        if (descr.initialValue().present()) {
            if (VERB >2) {
                std::cerr << "Initial direction for " << descr.name() << " is " << descr.initialValue().get() << std::endl;
            }
            init_values.push_back(descr.initialValue().get());
        }
    }
    
    FeatureNames& feNames = DataSet::getFeatureNames();
    feNames.setFromXml(featureDescriptions, DataSet::getCalcDoc());
    targetSet->initFeatureTables(feNames.getNumFeatures(), DataSet::getCalcDoc());
    decoySet->initFeatureTables(feNames.getNumFeatures(), DataSet::getCalcDoc());

    // import info from xml: read Fragment Spectrum Scans
    for (doc = p.next(); doc.get()!= 0 && 
          XMLString::equals(fragSpectrumScanStr, doc->getDocumentElement()->getTagName()); doc = p.next()) 
    {
      percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
      for (const auto &psm : fragSpectrumScan.peptideSpectrumMatch()) {
        if (psm.isDecoy()) {
          decoySet->readPsm(psm,fragSpectrumScan.scanNumber());
        } else {
          targetSet->readPsm(psm,fragSpectrumScan.scanNumber());
        }
      }
    }

    // import info from xml: read database proteins
    // only read them if they are present and the option of using mayusfdr is activated
    unsigned readProteins = 0;
    for (doc = p.next(); doc.get()!= 0 
        && hasProteins && ProteinProbEstimator::getCalcProteinLevelProb() /*&& Caller::protEstimator->getMayuFdr()*/
        && XMLString::equals(proteinStr, doc->getDocumentElement()->getTagName()); doc = p.next()) 
    {
      std::unique_ptr< ::percolatorInNs::protein > protein( new ::percolatorInNs::protein(*doc->getDocumentElement()));
      protEstimator->addProteinDb(*protein);
      ++readProteins;
    }
    
    /*if(ProteinProbEstimator::getCalcProteinLevelProb() && Caller::protEstimator->getMayuFdr() && readProteins <= 0)
    {
std::cerr << "Warning : options -Q and -A are activated but the number of proteins found in the input file is zero.\n\
	       Did you run converters with the flag -F ?\n" << std::endl;
Caller::protEstimator->setMayusFDR(false);
    }*/
    
    //maybe better to do :
    //SanityCheck::addDefaultWeights(init_values);
    pCheck = SanityCheck::initialize(otherCall);
    assert(pCheck);
    pCheck->addDefaultWeights(init_values);
    pCheck->checkAndSetDefaultDir();
    xmlInStream.close();
  } catch (const xml_schema::exception& e) {
    std::cerr << e << endl;
    return 0;
  } catch (const std::ios_base::failure&) {
    std::cerr << "ERROR: unable to open or read" << std::endl;
    return 0;
  } catch (const xercesc::DOMException& e) {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "ERROR: catched xercesc::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);
    return 0;
  }
  
  xercesc::XMLPlatformUtils::Terminate();

  setHandler.push_back_dataset(targetSet);
  setHandler.push_back_dataset(decoySet);
  return 1;
#else //XML_SUPPORT
  std::cerr << "ERROR: Compiler flag XML_SUPPORT was off, you cannot use the -k flag for pin-format input files" << std::endl;
  return 0;
#endif //XML_SUPPORT
}

/** 
 * Subroutine of @see Caller::writeXML() for PSM output
 */
void XMLInterface::writeXML_PSMs(Scores & fullset) {
  ofstream os;
  xmlOutputFN_PSMs = xmlOutputFN;
  xmlOutputFN_PSMs.append("writeXML_PSMs");
  os.open(xmlOutputFN_PSMs.c_str(), ios::out);

  os << "  <psms>" << endl;
  for (vector<ScoreHolder>::iterator psm = fullset.begin();
      psm != fullset.end(); ++psm) {
      os << *psm;
  }
  os << "  </psms>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see Caller::writeXML() for peptide output
 */
void XMLInterface::writeXML_Peptides(Scores & fullset) {
  reportUniquePeptides = true;
  ofstream os;
  xmlOutputFN_Peptides = xmlOutputFN;
  xmlOutputFN_Peptides.append("writeXML_Peptides");
  os.open(xmlOutputFN_Peptides.c_str(), ios::out);
  // append PEPTIDEs
  os << "  <peptides>" << endl;
  for (vector<ScoreHolder>::iterator psm = fullset.begin(); psm
  != fullset.end(); ++psm) {
    os << (ScoreHolderPeptide)*psm;
  }
  os << "  </peptides>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see Caller::writeXML() for protein output
 */
void XMLInterface::writeXML_Proteins(ProteinProbEstimator * protEstimator) {
  xmlOutputFN_Proteins = xmlOutputFN;
  xmlOutputFN_Proteins.append("writeXML_Proteins");
  protEstimator->writeOutputToXML(xmlOutputFN_Proteins, Scores::isOutXmlDecoys());
}

/** 
 * Writes the output of percolator to an pout XML file
 */
void XMLInterface::writeXML(Scores & fullset, ProteinProbEstimator * protEstimator, std::string call){
  ofstream os;
  const string space = PERCOLATOR_OUT_NAMESPACE;
  const string schema = space +
      " https://github.com/percolator/percolator/raw/pout-" + POUT_VERSION_MAJOR +
      "-" + POUT_VERSION_MINOR + "/src/xml/percolator_out.xsd";
  os.open(xmlOutputFN.data(), ios::out | ios::binary);
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  os << "<percolator_output "
      << endl << "xmlns=\""<< space << "\" "
      << endl << "xmlns:p=\""<< space << "\" "
      << endl << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
      << endl << "xsi:schemaLocation=\""<< schema <<"\" "
      << endl << "p:majorVersion=\"" << VERSION_MAJOR << "\" p:minorVersion=\""
      << VERSION_MINOR << "\" p:percolator_version=\"Percolator version "
      << VERSION << "\">\n"<< endl;
  os << "  <process_info>" << endl;
  os << "    <command_line>" << call << "</command_line>" << endl;

  os << "    <other_command_line>" << otherCall << "</other_command_line>\n";
  os << "    <pi_0_psms>" << pi_0_psms << "</pi_0_psms>" << endl;
  if(reportUniquePeptides)
    os << "    <pi_0_peptides>" << pi_0_peptides << "</pi_0_peptides>" << endl;
  if(ProteinProbEstimator::getCalcProteinLevelProb()) {  
    if(protEstimator->getUsePi0())
      os << "    <pi_0_proteins>" << protEstimator->getPi0() << "</pi_0_proteins>" << endl;
    /*if(protEstimator->getMayuFdr())
      os << "    <fdr_proteins>" << protEstimator->getFDR() << "</fdr_proteins>" << endl;*/
    os << "    <alpha>" << protEstimator->getAlpha() <<"</alpha>" << endl;
    os << "    <beta>"  << protEstimator->getBeta() <<"</beta>" << endl;
    os << "    <gamma>" << protEstimator->getGamma() <<"</gamma>" << endl;
  }
  os << "    <psms_qlevel>" <<  numberQpsms <<"</psms_qlevel>" << endl;
  if(reportUniquePeptides)
    os << "    <peptides_qlevel>" << fullset.getQvaluesBelowLevel(0.01) << "</peptides_qlevel>" << endl;
  if(ProteinProbEstimator::getCalcProteinLevelProb())
    os << "    <proteins_qlevel>" << protEstimator->getQvaluesBelowLevel(0.01) << "</proteins_qlevel>" << endl;  
  if (DataSet::getCalcDoc()) {
    os << "    <average_delta_mass>" << fullset.getDOC().getAvgDeltaMass()
                   << "</average_delta_mass>" << endl;
    os << "    <average_pi>" << fullset.getDOC().getAvgPI()
                   << "</average_pi>" << endl;
  }
  os << "  </process_info>" << endl << endl;

  // apppend PSMs
  ifstream ifs_psms(xmlOutputFN_PSMs.data(), ios::in | ios::binary);
  os << ifs_psms.rdbuf();
  ifs_psms.close();
  remove(xmlOutputFN_PSMs.c_str());
  // append Peptides
  if(reportUniquePeptides){
    ifstream ifs_peptides(xmlOutputFN_Peptides.data(), ios::in | ios::binary);
    os << ifs_peptides.rdbuf();
    ifs_peptides.close();
    remove(xmlOutputFN_Peptides.c_str());
  }
  // append Proteins
  if(ProteinProbEstimator::getCalcProteinLevelProb()){
    ifstream ifs_proteins(xmlOutputFN_Proteins.data(), ios::in | ios::binary);
    os << ifs_proteins.rdbuf();
    ifs_proteins.close();
    remove(xmlOutputFN_Proteins.c_str());
  }

  os << "</percolator_output>" << endl;
  os.close();
}

