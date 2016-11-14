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

#ifdef XML_SUPPORT

/** some constant strings to be used to compare xml strings **/

//databases
static const XMLCh databasesStr[] = {
  xercesc::chLatin_d, xercesc::chLatin_a, xercesc::chLatin_t, xercesc::chLatin_a, 
  xercesc::chLatin_b, xercesc::chLatin_a, xercesc::chLatin_s, xercesc::chLatin_e, 
  xercesc::chLatin_s, xercesc::chNull };
      
//calibration
static const XMLCh calibrationStr[] = { 
  xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_l, xercesc::chLatin_i, 
  xercesc::chLatin_b, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_t, 
  xercesc::chLatin_i, xercesc::chLatin_o, xercesc::chLatin_n, xercesc::chNull };

//proteins
static const XMLCh proteinsStr[] = { 
  xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t, 
  xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chLatin_s, 
  xercesc::chNull };

//protein      
static const XMLCh proteinStr[] = { 
  xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t, 
  xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chNull };
      
//fragSpectrumScan 
static const XMLCh fragSpectrumScanStr[] = { 
  xercesc::chLatin_f, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_g, 
  xercesc::chLatin_S, xercesc::chLatin_p, xercesc::chLatin_e, xercesc::chLatin_c, 
  xercesc::chLatin_t, xercesc::chLatin_r, xercesc::chLatin_u, xercesc::chLatin_m,
  xercesc::chLatin_S, xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_n, 
  xercesc::chNull };  
      
#endif //XML_SUPPORT

XMLInterface::XMLInterface(const std::string& outputFN, 
    bool schemaValidation, bool printDecoys, bool printExpMass) : 
  xmlOutputFN_(outputFN), schemaValidation_(schemaValidation), 
  otherCall_(""), reportUniquePeptides_(false), printDecoys_(printDecoys), 
  printExpMass_(printExpMass) {}

XMLInterface::~XMLInterface() {
  // clean up temporary files if an exception occurred during the writing
  remove(xmlOutputFN_PSMs.c_str());
  remove(xmlOutputFN_Peptides.c_str());
  remove(xmlOutputFN_Proteins.c_str());
}

int XMLInterface::readPin(istream& dataStream, const std::string& xmlInputFN,
    SetHandler& setHandler, SanityCheck*& pCheck, 
    ProteinProbEstimator* protEstimator) {
  std::vector<double> noWeights;
  Scores noScores(true);
  return readAndScorePin(dataStream, noWeights, noScores, xmlInputFN, setHandler, 
                  pCheck, protEstimator);
}
  
int XMLInterface::readAndScorePin(istream& dataStream, std::vector<double>& rawWeights, 
    Scores& allScores, const std::string& xmlInputFN,
    SetHandler& setHandler, SanityCheck*& pCheck, 
    ProteinProbEstimator* protEstimator) {    
#ifdef XML_SUPPORT
  xercesc::XMLPlatformUtils::Initialize();
  try {
    using namespace xercesc;

    string schemaDefinition = Globals::getInstance()->getXMLDir()+PIN_SCHEMA_LOCATION+string("percolator_in.xsd");
    parser p;
    xml_schema::dom::auto_ptr<DOMDocument> doc(p.start(
        dataStream, xmlInputFN.c_str(), schemaValidation_,
        schemaDefinition, PIN_VERSION_MAJOR, PIN_VERSION_MINOR));

    doc = p.next();
    // read enzyme element
    // the enzyme element is a subelement but CodeSynthesis Xsd does not
    // generate a class for it. (I am trying to find a command line option
    // that overrides this decision). As for now special treatment is needed
    char* value = XMLString::transcode(doc->getDocumentElement()->getTextContent());
    
    if (VERB > 1) std::cerr << "enzyme=" << value << std::endl;
    
    Enzyme::setEnzyme(value);
    XMLString::release(&value);
    doc = p.next();
    
    //checking if database is present to jump it
    bool hasProteins = false;
    if (XMLString::equals(databasesStr, doc->getDocumentElement()->getTagName())) {
      //NOTE I dont really need this info, do I? good to have it though
      // std::unique_ptr< ::percolatorInNs::databases > databases( new ::percolatorInNs::databases(*doc->getDocumentElement()));
      doc = p.next();
      hasProteins = true;
    }
    
    // read process_info element
    percolatorInNs::process_info processInfo(*doc->getDocumentElement());
    otherCall_ = processInfo.command_line();
    doc = p.next();

    if (XMLString::equals(calibrationStr,doc->getDocumentElement()->getTagName())) {
      //NOTE the calibration should define the initial direction
      //percolatorInNs::calibration calibration(*doc->getDocumentElement());
      doc = p.next();
    };
    
    // read feature names and initial values that are present in feature descriptions
    FeatureNames& featureNames = DataSet::getFeatureNames();
    percolatorInNs::featureDescriptions featureDescriptions(*doc->getDocumentElement());
    percolatorInNs::featureDescriptions::featureDescription_const_iterator featureIt;
    featureIt = featureDescriptions.featureDescription().begin();
    for ( ; featureIt != featureDescriptions.featureDescription().end(); ++featureIt) {
      featureNames.insertFeature(featureIt->name());
    }
    featureNames.initFeatures(DataSet::getCalcDoc());
    
    std::vector<double> init_values(FeatureNames::getNumFeatures());
    setHandler.getFeaturePool().createPool(FeatureNames::getNumFeatures());
    
    bool hasDefaultValues = false;
    unsigned int i = 0;
    featureIt = featureDescriptions.featureDescription().begin();
    for ( ; featureIt != featureDescriptions.featureDescription().end(); ++featureIt) {    
      if (featureIt->initialValue().present()) {
        if (featureIt->initialValue().get() != 0.0) 
          hasDefaultValues = true;
        if (VERB > 2) {
          std::cerr << "Initial direction for " << featureIt->name() << " is " << 
                       featureIt->initialValue().get() << std::endl;
        }
        init_values[i] = featureIt->initialValue().get();
      }
      ++i;
    }
    
    bool readProteins = true;
    if (rawWeights.size() == 0) {
      DataSet* targetSet = new DataSet();
      assert(targetSet);
      targetSet->setLabel(1);
      DataSet* decoySet = new DataSet();
      assert(decoySet);
      decoySet->setLabel(-1);    
      
      // detect if the input came from separate target and decoy searches or 
      // from a concatenated search by looking for scan+expMass combinations
      // that have both at least one target and decoy PSM
      bool concatenatedSearch = true;
    
      // read Fragment Spectrum Scans
      if (setHandler.getMaxPSMs() > 0u) { // subset training option is turned on
        readProteins = false;
        std::priority_queue<PSMDescriptionPriority> subsetPSMs;
        std::map<ScanId, std::pair<size_t, bool> > scanIdLookUp; // ScanId -> (priority, isDecoy)
        unsigned int upperLimit = UINT_MAX;
        for (doc = p.next(); 
             doc.get()!= 0 && XMLString::equals(fragSpectrumScanStr, 
                 doc->getDocumentElement()->getTagName()); 
             doc = p.next()) {
          percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
          percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_const_iterator psmIt;
          psmIt = fragSpectrumScan.peptideSpectrumMatch().begin();
          for ( ; psmIt != fragSpectrumScan.peptideSpectrumMatch().end(); ++psmIt) {
            ScanId scanId = getScanId(*psmIt, fragSpectrumScan.scanNumber());
            size_t randIdx;
            if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
              if (concatenatedSearch && psmIt->isDecoy() != scanIdLookUp[scanId].second) {
                concatenatedSearch = false;
              }
              randIdx = scanIdLookUp[scanId].first;
            } else {
              randIdx = PseudoRandom::lcg_rand();
              scanIdLookUp[scanId].first = randIdx;
              scanIdLookUp[scanId].second = psmIt->isDecoy();
            }
            
            if (subsetPSMs.size() < setHandler.getMaxPSMs() || randIdx < upperLimit) {
              PSMDescriptionPriority psmPriority;
              psmPriority.psm = readPsm(*psmIt, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());
              psmPriority.label = (psmIt->isDecoy() ? -1 : 1);
              psmPriority.priority = randIdx;
              subsetPSMs.push(psmPriority);
              if (subsetPSMs.size() > setHandler.getMaxPSMs()) {
                PSMDescriptionPriority del = subsetPSMs.top();
                upperLimit = del.priority;
                setHandler.getFeaturePool().deallocate(del.psm->features);
                PSMDescription::deletePtr(del.psm);
                subsetPSMs.pop();
              }
            }
          }
        }
        
        setHandler.addQueueToSets(subsetPSMs, targetSet, decoySet);
      } else {
        std::map<ScanId, bool> scanIdLookUp; // ScanId -> isDecoy
        for (doc = p.next(); 
             doc.get()!= 0 && XMLString::equals(fragSpectrumScanStr, 
                 doc->getDocumentElement()->getTagName()); 
             doc = p.next()) {
          percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
          percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_const_iterator psmIt;
          psmIt = fragSpectrumScan.peptideSpectrumMatch().begin();
          for ( ; psmIt != fragSpectrumScan.peptideSpectrumMatch().end(); ++psmIt) {
            ScanId scanId = getScanId(*psmIt, fragSpectrumScan.scanNumber());
            if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
              if (concatenatedSearch && psmIt->isDecoy() != scanIdLookUp[scanId]) {
                concatenatedSearch = false;
              }
            } else {
              scanIdLookUp[scanId] = psmIt->isDecoy();
            }
            
            PSMDescription* psm = readPsm(*psmIt, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());
            if (psmIt->isDecoy()) {
              decoySet->registerPsm(psm);
            } else {
              targetSet->registerPsm(psm);
            }
          }
        }
      }
    
      setHandler.push_back_dataset(targetSet);
      setHandler.push_back_dataset(decoySet);
      
      pCheck = SanityCheck::initialize(otherCall_);
      assert(pCheck);
      pCheck->checkAndSetDefaultDir();
      if (hasDefaultValues) pCheck->addDefaultWeights(init_values);
      pCheck->setConcatenatedSearch(concatenatedSearch);
      
    } else { // read PSMs and immediately score them (second step of subset training option)
      const unsigned int numFeatures = FeatureNames::getNumFeatures();
      for (doc = p.next(); 
           doc.get()!= 0 && XMLString::equals(fragSpectrumScanStr, 
               doc->getDocumentElement()->getTagName()); 
           doc = p.next()) {
        percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
        percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_const_iterator psmIt;
        psmIt = fragSpectrumScan.peptideSpectrumMatch().begin();
        for ( ; psmIt != fragSpectrumScan.peptideSpectrumMatch().end(); ++psmIt) {
          ScoreHolder sh;
          sh.label = (psmIt->isDecoy() ? -1 : 1);
          sh.pPSM = readPsm(*psmIt, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());
          
          allScores.scoreAndAddPSM(sh, rawWeights, setHandler.getFeaturePool());
        }
      }
    }
    
    if (readProteins) {
      // read database proteins
      // only read them if they are present and the option of using mayusfdr is activated
      unsigned numProteins = 0;
      if (hasProteins && ProteinProbEstimator::getCalcProteinLevelProb() /*&& Caller::protEstimator->getMayuFdr()*/) {
        assert(protEstimator); // should be initialized if -A or -f flag was used (which also sets calcProteinLevelProb)
        for (doc = p.next(); doc.get()!= 0
            && XMLString::equals(proteinStr, doc->getDocumentElement()->getTagName()); doc = p.next()) 
        {
          ::percolatorInNs::protein protein(*doc->getDocumentElement());
          protEstimator->addProteinDb(protein.isDecoy(), protein.name(), protein.sequence(), protein.length());
          ++numProteins;
        }
      }
      /*
      if(ProteinProbEstimator::getCalcProteinLevelProb() && protEstimator->getMayuFdr() && numProteins <= 0) {
        std::cerr << "Warning : options -Q and -A are activated but the number of proteins found in the input file is zero.\n\
	         Did you run converters with the flag -F ?\n" << std::endl;
        Caller::protEstimator->setMayusFDR(false);
      }
      */
    }
    
  } catch (const xml_schema::exception& e) {
    std::cerr << "ERROR: xml schema error " << e << endl;
    return 0;
  } catch (const std::ios_base::failure&) {
    std::cerr << "ERROR: unable to open or read" << std::endl;
    return 0;
  } catch (const xercesc::DOMException& e) {
    char* tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "ERROR: caught xercesc::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);
    return 0;
  }
  
  xercesc::XMLPlatformUtils::Terminate();
  return 1;
#else //XML_SUPPORT
  std::cerr << "ERROR: Compiler flag XML_SUPPORT was off, you cannot use the -k flag for pin-format input files" << std::endl;
  return 0;
#endif //XML_SUPPORT
}


#ifdef XML_SUPPORT
// Convert a peptide with or without modifications into a string
std::string XMLInterface::decoratePeptide(const ::percolatorInNs::peptideType& peptide) {
  std::list<std::pair<int,std::string> > mods;
  std::string peptideSeq = peptide.peptideSequence();
  percolatorInNs::peptideType::modification_const_iterator modIt;
  modIt = peptide.modification().begin();
  for ( ; modIt != peptide.modification().end(); ++modIt) {
    std::stringstream ss;
    if (modIt->uniMod().present()) {
      ss << "[UNIMOD:" << modIt->uniMod().get().accession() << "]";
      mods.push_back(std::pair<int,std::string>(modIt->location(),ss.str()));
    }
    if (modIt->freeMod().present()) {
      ss << "[" << modIt->freeMod().get().moniker() << "]";
      mods.push_back(std::pair<int,std::string>(modIt->location(),ss.str()));
    }
  }
  mods.sort(greater<std::pair<int,std::string> >());
  std::list<std::pair<int,std::string> >::const_iterator it;
  for(it=mods.begin();it!=mods.end();++it) {
    peptideSeq.insert(it->first,it->second);
  }
  return peptideSeq;
}

PSMDescription* XMLInterface::readPsm(
    const percolatorInNs::peptideSpectrumMatch& psm, unsigned scanNumber, 
    bool readProteins, FeatureMemoryPool& featurePool) {
  PSMDescription* myPsm;
  if (DataSet::getCalcDoc()) {
    myPsm = new PSMDescriptionDOC();
  } else {
    myPsm = new PSMDescription();
  }
  string mypept = decoratePeptide(psm.peptide());

  if (psm.occurence().size() <= 0) {
    ostringstream temp;
    temp << "Error: cannot add PSM " << psm.id() << " to the dataset.\n\
    The PSM does not contain protein occurences." << std::endl;
    throw MyException(temp.str());
  }
  
  percolatorInNs::peptideSpectrumMatch::occurence_const_iterator occIt;
  occIt = psm.occurence().begin();
  for ( ; occIt != psm.occurence().end(); ++occIt) {
    if (readProteins) myPsm->proteinIds.push_back( occIt->proteinId() );
    // adding n-term and c-term residues to peptide
    //NOTE the residues for the peptide in the PSMs are always the same for every protein
    myPsm->peptide = occIt->flankN() + "." + mypept + "." + occIt->flankC();
  }

  myPsm->setId(psm.id());
  myPsm->scan = scanNumber;
  myPsm->expMass = psm.experimentalMass();
  myPsm->calcMass = psm.calculatedMass();
  if ( psm.observedTime().present() ) {
    myPsm->setRetentionTime(psm.observedTime().get());
  }

  myPsm->features = featurePool.allocate();
  
  for (unsigned int i = 0; i < psm.features().feature().size(); ++i) {
    myPsm->features[i] = psm.features().feature()[i];
  }

  // myPsm.peptide = psmIter->peptide().peptideSequence();
  myPsm->setMassDiff(MassHandler::massDiff(psm.experimentalMass(), 
      psm.calculatedMass(), psm.chargeState()));
  return myPsm;
}

ScanId XMLInterface::getScanId(const percolatorInNs::peptideSpectrumMatch& psm, 
                               unsigned scanNumber) {
  ScanId scanId;
  scanId.first = scanNumber;
  scanId.second = psm.experimentalMass();
  return scanId;
}
#endif // XML_SUPPORT

/** 
 * Subroutine of @see XMLInterface::writeXML() for PSM output
 */
void XMLInterface::writeXML_PSMs(Scores& fullset) {
  pi0Psms_ = fullset.getPi0();
  numberQpsms_ = fullset.getQvaluesBelowLevel(0.01);
  
  ofstream os;
  xmlOutputFN_PSMs = xmlOutputFN_;
  xmlOutputFN_PSMs.append("writeXML_PSMs");
  os.open(xmlOutputFN_PSMs.c_str(), ios::out);
  
  os << "  <psms>" << endl;
  for (std::vector<ScoreHolder>::iterator psm = fullset.begin();
       psm != fullset.end(); ++psm) {
    psm->printPSM(os, printDecoys_, printExpMass_);
  }
  os << "  </psms>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see XMLInterface::writeXML() for peptide output
 */
void XMLInterface::writeXML_Peptides(Scores& fullset) {
  pi0Peptides_ = fullset.getPi0();
  reportUniquePeptides_ = true;
  
  ofstream os;
  xmlOutputFN_Peptides = xmlOutputFN_;
  xmlOutputFN_Peptides.append("writeXML_Peptides");
  os.open(xmlOutputFN_Peptides.c_str(), ios::out);
  // append PEPTIDEs
  os << "  <peptides>" << endl;
  for (vector<ScoreHolder>::iterator psm = fullset.begin(); 
       psm != fullset.end(); ++psm) {
    psm->printPeptide(os, printDecoys_, printExpMass_, fullset);
  }
  os << "  </peptides>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see XMLInterface::writeXML() for protein output
 */
void XMLInterface::writeXML_Proteins(ProteinProbEstimator * protEstimator) {
  xmlOutputFN_Proteins = xmlOutputFN_;
  xmlOutputFN_Proteins.append("writeXML_Proteins");
  protEstimator->writeOutputToXML(xmlOutputFN_Proteins, printDecoys_);
}

/** 
 * Writes the output of percolator to an pout XML file
 */
void XMLInterface::writeXML(Scores& fullset, ProteinProbEstimator* protEstimator, std::string call) {
  ofstream os;
  const string space = PERCOLATOR_OUT_NAMESPACE;
  const string schema = space +
      " https://github.com/percolator/percolator/raw/pout-" + POUT_VERSION_MAJOR +
      "-" + POUT_VERSION_MINOR + "/src/xml/percolator_out.xsd";
  os.open(xmlOutputFN_.data(), ios::out | ios::binary);
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

  os << "    <other_command_line>" << otherCall_ << "</other_command_line>\n";
  os << "    <pi_0_psms>" << pi0Psms_ << "</pi_0_psms>" << endl;
  if (reportUniquePeptides_)
    os << "    <pi_0_peptides>" << pi0Peptides_ << "</pi_0_peptides>" << endl;
  if (ProteinProbEstimator::getCalcProteinLevelProb()) {  
    if (protEstimator->getUsePi0())
      os << "    <pi_0_proteins>" << protEstimator->getPi0() << "</pi_0_proteins>" << endl;
    /*if(protEstimator->getMayuFdr())
      os << "    <fdr_proteins>" << protEstimator->getFDR() << "</fdr_proteins>" << endl;*/
    protEstimator->printParametersXML(os);
  }
  os << "    <psms_qlevel>" <<  numberQpsms_ <<"</psms_qlevel>" << endl;
  if (reportUniquePeptides_)
    os << "    <peptides_qlevel>" << fullset.getQvaluesBelowLevel(0.01) << "</peptides_qlevel>" << endl;
  if (ProteinProbEstimator::getCalcProteinLevelProb())
    os << "    <proteins_qlevel>" << protEstimator->getQvaluesBelowLevel(0.01) << "</proteins_qlevel>" << endl;  
  if (DataSet::getCalcDoc()) {
    os << "    <average_delta_mass>" << fullset.getDOC().getAvgDeltaMass()
                   << "</average_delta_mass>" << endl;
    os << "    <average_pi>" << fullset.getDOC().getAvgPI()
                   << "</average_pi>" << endl;
  }
  os << "  </process_info>" << endl << endl;

  // append PSMs
  ifstream ifs_psms(xmlOutputFN_PSMs.data(), ios::in | ios::binary);
  os << ifs_psms.rdbuf();
  ifs_psms.close();
  remove(xmlOutputFN_PSMs.c_str());
  // append Peptides
  if (reportUniquePeptides_){
    ifstream ifs_peptides(xmlOutputFN_Peptides.data(), ios::in | ios::binary);
    os << ifs_peptides.rdbuf();
    ifs_peptides.close();
    remove(xmlOutputFN_Peptides.c_str());
  }
  // append Proteins
  if (ProteinProbEstimator::getCalcProteinLevelProb()){
    ifstream ifs_proteins(xmlOutputFN_Proteins.data(), ios::in | ios::binary);
    os << ifs_proteins.rdbuf();
    ifs_proteins.close();
    remove(xmlOutputFN_Proteins.c_str());
  }
  os << "</percolator_output>" << endl;
  os.close();
}

