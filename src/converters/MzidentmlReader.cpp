#include <memory> // For std::unique_ptr
#include <xercesc/dom/DOM.hpp> // For xercesc::DOMDocument
#include <xsd/cxx/xml/dom/auto-ptr.hxx> // For xsd::cxx::xml::dom::deleter

#include "MzidentmlReader.h"

static const XMLCh sequenceCollectionStr[] = {
    chLatin_S, chLatin_e, chLatin_q, chLatin_u, chLatin_e, chLatin_n,
    chLatin_c, chLatin_e, chLatin_C, chLatin_o, chLatin_l, chLatin_l, chLatin_e,
    chLatin_c, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull};


static string schemaDefinition = Globals::getInstance()->getXMLDir(true) + 
           MZIDENTML_SCHEMA_LOCATION + string("mzIdentML1.1.0.xsd");
static string scheme_namespace = MZIDENTML_NAMESPACE;
static string schema_major = boost::lexical_cast<string>(MZIDENTML_VERSION_MAJOR);
static string schema_minor = boost::lexical_cast<string>(MZIDENTML_VERSION_MINOR);

static double proton_mass = 1.00727663;

MzidentmlReader::MzidentmlReader(ParseOptions po) : Reader(po) {}

MzidentmlReader::~MzidentmlReader() {}

void MzidentmlReader::cleanHashMaps() {
  peptideMapType::iterator iter;
  for (iter = peptideMap.begin(); iter != peptideMap.end(); ++iter) {
    if(iter->second) delete iter->second;
    iter->second = 0;
  }

  proteinMapType::iterator iter2;
  for (iter2 = proteinMap.begin(); iter2 != proteinMap.end(); ++iter2) {
    if(iter2->second) delete iter2->second;
    iter2->second = 0;
  }

  peptideEvidenceMapType::iterator iter3;
  for (iter3 = peptideEvidenceMap.begin(); iter3 != peptideEvidenceMap.end(); ++iter3) {
    if(iter3->second) delete iter3->second;
    iter3->second = 0;
  }
}

bool MzidentmlReader::checkIsMeta(const std::string &file) {
  //NOTE assuming the file has been tested before
  bool isMeta;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  std::string line;
  getline(fileIn, line);
  fileIn.close();
  //NOTE this is not the best way to check if it is meta for mzident
  if (line.find("<?xml") != std::string::npos) {
    isMeta = false;
  } else {
    isMeta = true;
  }
  return isMeta;
}

void MzidentmlReader::getMaxMinCharge(const std::string &fn, bool isDecoy) {
  ifstream ifs;
  ifs.exceptions(ifstream::badbit | ifstream::failbit);
  try {
    ifs.open(fn.c_str());
    parser p;
    bool validateSchema = true;

    std::unique_ptr<xercesc::DOMDocument, xsd::cxx::xml::dom::deleter<xercesc::DOMDocument>> doc(nullptr, xsd::cxx::xml::dom::deleter<xercesc::DOMDocument>());

    // Ensure parser::start returns the correct type or wrap the result with the correct deleter
    std::unique_ptr<xercesc::DOMDocument> tempDoc = p.start(ifs, fn.c_str(), validateSchema, schemaDefinition, schema_major, schema_minor, scheme_namespace);
    doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter

    // MT: This seems to be a bit slow for doing nothing
    for (std::unique_ptr<xercesc::DOMDocument> tempDoc = p.next(); tempDoc.get() != nullptr
        && !XMLString::equals(spectrumIdentificationResultStr, tempDoc->getDocumentElement()->getTagName()); tempDoc = p.next()) {
      doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter
    }
    
    // For each SpectrumIdentificationResult
    int itemCount = 1;
    for (auto tempDoc = p.next(); tempDoc.get() != nullptr && XMLString::equals(spectrumIdentificationResultStr,
            tempDoc->getDocumentElement()->getTagName()); tempDoc = p.next()) {
      doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter
      ::mzIdentML_ns::SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement());
      // For each SpectrumIdentificationItem
      BOOST_FOREACH(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem()) {
        minCharge = (std::min)(item.chargeState(), minCharge);
        maxCharge = (std::max)(item.chargeState(), maxCharge);
        searchEngineSpecificParsing(item, itemCount);  // Virtual function that potentially checks the features
        ++itemCount;
      }
    }
  } catch (ifstream::failure e) {
    std::cerr << "Exception opening/reading file:" << fn << std::endl;
  } catch (const xercesc::DOMException& e) {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "Error parsing file: " << fn << "\n  xercesc::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);
  } catch (const xml_schema::exception& e) {
    std::cerr << "Error parsing file: " << fn << "\n  XML schema exception in getMaxMinCharge: " << e << std::endl;
  } catch (MyException e) {
	  std::cerr << "Error reading file: " << fn << "\n  " << e.what() << std::endl;
	  exit(1);
  } catch (std::exception e) {
    std::cerr << "Error reading file: " << fn << "\n  Unknown exception in getMaxMinCharge: " << e.what() << std::endl;
  }
  ifs.close();
  return;
}


void MzidentmlReader::searchEngineSpecificParsing(
    const ::mzIdentML_ns::SpectrumIdentificationItemType & item, 
    int itemCount) {
	// itemCount is which order the item is in the mzid-file, read through userParam elements
	return;  // Empty function
}

void MzidentmlReader::read(const std::string &fn, bool isDecoy, 
    boost::shared_ptr<FragSpectrumScanDatabase> database) {
  namespace xml = xsd::cxx::xml;
  scanNumberMapType scanNumberMap;
  ifstream ifs;
  try
  {
    ifs.exceptions(ifstream::badbit | ifstream::failbit);
    ifs.open(fn.c_str());
    parser p;
    
    std::unique_ptr<xercesc::DOMDocument, xsd::cxx::xml::dom::deleter<xercesc::DOMDocument>> doc;

    // Ensure parser::start returns the correct type or wrap the result with the correct deleter
    auto tempDoc = p.start(ifs, fn.c_str(), true, schemaDefinition, schema_major, schema_minor, scheme_namespace);
    doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter

    //NOTE wouldnt be  better to use the get tag by Name to jump SequenceCollenction directly?
    while (doc.get() != nullptr && !XMLString::equals(sequenceCollectionStr,
            doc->getDocumentElement()->getTagName())) {
      auto tempDoc = p.next();
      doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter
    }

    assert(doc.get());
    mzIdentML_ns::SequenceCollectionType sequenceCollection(*doc->getDocumentElement());

    peptideMap.clear();
    proteinMap.clear();
    peptideEvidenceMap.clear();

    //NOTE probably I can get rid of these hash tables with a proper access to elements by tag and id

    BOOST_FOREACH(const mzIdentML_ns::SequenceCollectionType::Peptide_type &peptide, sequenceCollection.Peptide()) {
      //PEPTIDE
      mzIdentML_ns::SequenceCollectionType::Peptide_type *pept =
              new mzIdentML_ns::SequenceCollectionType::Peptide_type(peptide);
      peptideMap.insert(std::make_pair(peptide.id(), pept));
    }

    BOOST_FOREACH(const mzIdentML_ns::SequenceCollectionType::DBSequence_type &protein, sequenceCollection.DBSequence()) {
      //PROTEIN
      mzIdentML_ns::SequenceCollectionType::DBSequence_type *prot =
              new mzIdentML_ns::SequenceCollectionType::DBSequence_type(protein);
      proteinMap.insert(std::make_pair(protein.id(), prot));
    }

    BOOST_FOREACH(const ::mzIdentML_ns::PeptideEvidenceType &peptideE, sequenceCollection.PeptideEvidence()) {
      //PEPTIDE EVIDENCE
      ::mzIdentML_ns::PeptideEvidenceType *peptE = new mzIdentML_ns::PeptideEvidenceType(peptideE);
      peptideEvidenceMap.insert(std::make_pair(peptideE.id(), peptE));
    }

    for (auto tempDoc = p.next(); tempDoc.get() != nullptr && !XMLString::equals(spectrumIdentificationResultStr,
            tempDoc->getDocumentElement()->getTagName()); tempDoc = p.next()) {
      doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter
    }

    unsigned scanNumber = 0;
    bool useRankedScanNumbers = false;  	/* True scan numbers are used,
    										if they can't be found, use ranked scan numbers from 1 and up. */
    for (auto tempDoc = p.next(); tempDoc.get() != nullptr && XMLString::equals(spectrumIdentificationResultStr,
            tempDoc->getDocumentElement()->getTagName()); tempDoc = p.next()) {
      doc.reset(tempDoc.release()); // Transfer ownership and apply the custom deleter
      ::mzIdentML_ns::SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement());
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      int numberHitsSpectra = 0;

      //Find scan number from the cvParam element in spetrumIdentificationResults
      if(!useRankedScanNumbers) {
    	  bool foundScanNumber = false;  // Indicates whether a proper scan number was found
    	  BOOST_FOREACH(const ::mzIdentML_ns::CVParamType & cv, specIdResult.cvParam()) {
    		  std::string param_name(cv.name().c_str());
    		  std::string expected_name("scan number(s)");
    		  if (param_name == expected_name) {
    			  scanNumber = boost::lexical_cast<unsigned>(cv.value().get().c_str());
    			  foundScanNumber = true;
    		  }
    	  }
    	  if(!foundScanNumber || scanNumber == 0) {
    		  std::cerr << "No scan number was found for a PSM (or it equaled 0), scans are ranked from 1 and up" << std::endl;
    		  useRankedScanNumbers = true;
    	  }
      }
      // If no scan numbers were found, or a scan of 0 was found, just rank them
      if(useRankedScanNumbers) {
    	  ++scanNumber;
      }

      BOOST_FOREACH(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem()) {
	      if(++numberHitsSpectra <= po.hitsPerSpectrum) {
	        assert(item.experimentalMassToCharge());
          int charge = item.chargeState();
	        ::percolatorInNs::fragSpectrumScan::experimentalMass_type experimentalMass = item.experimentalMassToCharge()*charge - proton_mass*charge;
	        createPSM(item, experimentalMass, isDecoy, scanNumber, database, fn);
	      }
      }
    }

    cleanHashMaps();
    ifs.close();
  }
  catch (const xercesc::DOMException& e)
  {
    cleanHashMaps();
    ifs.close();
    char * tmpStr = XMLString::transcode(e.getMessage());
    ostringstream temp;
    temp << "Error : xercesc::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);
    throw MyException(temp.str());
  }

  return;
}