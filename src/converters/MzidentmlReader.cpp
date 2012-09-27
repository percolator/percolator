#include "MzidentmlReader.h"

static const XMLCh sequenceCollectionStr[] = {chLatin_S, chLatin_e, chLatin_q, chLatin_u, chLatin_e, chLatin_n,
  chLatin_c, chLatin_e, chLatin_C, chLatin_o, chLatin_l, chLatin_l, chLatin_e,
  chLatin_c, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull};

//Hash maps between feature names indicated in the mzident-files
const std::map<string, int> MzidentmlReader::sequestFeatures =
        boost::assign::map_list_of("sequest:PeptideRankSp", 0)
                                  ("sequest:deltacn", 1)
                                  ("sequest:xcorr", 2)
                                  ("sequest:PeptideSp", 3)
                                  ("sequest:matched ions", 4)
                                  ("sequest:total ions", 5)
                                  ("sequest:PeptideIdnumber", 6)
                                  ("sequest:PeptideNumber", 7);

const std::map<string, int> MzidentmlReader::msgfplusFeatures =
        boost::assign::map_list_of("MS-GF:RawScore", 0)
                                  ("MS-GF:DeNovoScore", 1)
                                  ("MS-GF:SpecEValue", 2)
                                  ("MS-GF:EValue", 3)
                                  //The below features are on user specified element userParam
                                  ("IsotopeError", 4) 
                                  ("ExplainedIonCurrentRatio", 5)
                                  ("NTermIonCurrentRatio", 6)
                                  ("CTermIonCurrentRatio", 7)
                                  ("MS2IonCurrent", 8);

MzidentmlReader::MzidentmlReader(ParseOptions *po) : Reader(po) {

}

MzidentmlReader::~MzidentmlReader() {



}

//enum MzidentmlReader::inputFormat_t {sequest, msgfplus};  //These formats are accepted

void MzidentmlReader::cleanHashMaps() {
  peptideMapType::iterator iter;
  for (iter = peptideMap.begin(); iter != peptideMap.end(); ++iter) {
    delete iter->second;
    iter->second = 0;
  }

  proteinMapType::iterator iter2;
  for (iter2 = proteinMap.begin(); iter2 != proteinMap.end(); ++iter2) {
    delete iter2->second;
    iter2->second = 0;
  }

  peptideEvidenceMapType::iterator iter3;
  for (iter3 = peptideEvidenceMap.begin(); iter3 != peptideEvidenceMap.end(); ++iter3) {
    delete iter3->second;
    iter3->second = 0;
  }
}

bool MzidentmlReader::checkValidity(const std::string &file) {
  bool isvalid = true;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  if (!fileIn) {
    std::cerr << "Could not open file " << file << std::endl;
    exit(-1);
  }
  std::string line;
  if (!getline(fileIn, line)) {
    std::cerr << "Could not read file " << file << std::endl;
    exit(-1);
  }
  if (line.find("<?xml") == std::string::npos) {
    std::cerr << "ERROR : the input file is not xml format " << file << std::endl;
    exit(-1);
  } 
  else //Test whether Sequest or MS-GF+ format
  {
    std::string line2, line3;
    getline(fileIn, line2);
    getline(fileIn, line3);

    if ((line2[1] != '!' && line2.find("SEQUEST") != std::string::npos && line2.find("MzIdentML") != std::string::npos)
            || (line3[1] != '!' && line3.find("SEQUEST") != std::string::npos && line3.find("MzIdentML") != std::string::npos)) {
      std::cerr << "MzIdentML - Sequest format " << std::endl;
      inputFormat = sequest;

    } else if ((line2[1] != '!' && line2.find("MS-GF+") != std::string::npos && line2.find("MzIdentML") != std::string::npos)
            || (line3[1] != '!' && line3.find("MS-GF+") != std::string::npos && line3.find("MzIdentML") != std::string::npos)) {
      std::cerr << "MzIdentML - MSGF+ format" << std::endl;
      inputFormat = msgfplus;

    } else {
      std::cerr << "ERROR : the input file is not MzIdentML - Sequest or MSGF+ format " << file << std::endl;
      exit(-1);
    }

  }
  fileIn.close();
  return isvalid;
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

void MzidentmlReader::addFeatureDescriptions(bool doEnzyme) {
  if (inputFormat == sequest) 
  {
    push_backFeatureDescription("lnrSp");
    push_backFeatureDescription("deltCn");
    push_backFeatureDescription("Xcorr");
    push_backFeatureDescription("Sp");
    push_backFeatureDescription("IonFrac");
  } 
  else if (inputFormat == msgfplus) 
  {
    push_backFeatureDescription("RawScore");
    push_backFeatureDescription("DeNovoScore");
    push_backFeatureDescription("SpecEValue");
    push_backFeatureDescription("EValue");
    //The below are from element userParam
    push_backFeatureDescription("IsotopeError");
    push_backFeatureDescription("ExplainedIonCurrentRatio");
    push_backFeatureDescription("NTermIonCurrentRatio");
    push_backFeatureDescription("CTermIonCurrentRatio");
    push_backFeatureDescription("MS2IonCurrent");
  }
  push_backFeatureDescription("Mass");
  push_backFeatureDescription("PepLen");
  for (int charge = minCharge; charge <= maxCharge; ++charge) {
    std::ostringstream cname;
    cname << "Charge" << charge;
    push_backFeatureDescription(cname.str().c_str());

  }
  if (doEnzyme) {
    push_backFeatureDescription("enzN");
    push_backFeatureDescription("enzC");
    push_backFeatureDescription("enzInt");
  }

  push_backFeatureDescription("dM");
  push_backFeatureDescription("absdM");

  if (po->calcPTMs) {
    push_backFeatureDescription("ptm");
  }
  if (po->pngasef) {
    push_backFeatureDescription("PNGaseF");
  }

  if (po->calcAAFrequencies) {
    for (std::string::const_iterator it = aaAlphabet.begin(); it != aaAlphabet.end(); it++) {
      std::string temp = boost::lexical_cast<std::string > (*it) + "-Freq";
      push_backFeatureDescription(temp.c_str());
    }
  }

}

void MzidentmlReader::getMaxMinCharge(const std::string &fn, bool isDecoy) {

  ifstream ifs;
  ifs.exceptions(ifstream::badbit | ifstream::failbit);
  try {
    ifs.open(fn.c_str());
    parser p;
    string schemaDefinition = MZIDENTML_SCHEMA_LOCATION + string("mzIdentML1.1.0.xsd");
    string scheme_namespace = MZIDENTML_NAMESPACE;
    string schema_major = "";
    string schema_minor = "";
    bool schemaVal = true;
    xml_schema::dom::auto_ptr<DOMDocument> doc
            (p.start(ifs, fn.c_str(), schemaVal, schemaDefinition, schema_major, schema_minor, scheme_namespace));

    for (doc = p.next(); doc.get() != 0
            && !XMLString::equals(spectrumIdentificationResultStr, doc->getDocumentElement()->getTagName()); doc = p.next()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }

    for (; doc.get() != 0 && XMLString::equals(spectrumIdentificationResultStr,
            doc->getDocumentElement()->getTagName()); doc = p.next()) {
      ::mzIdentML_ns::SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement());

      BOOST_FOREACH(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem()) {
        minCharge = std::min(item.chargeState(), minCharge);
        maxCharge = std::max(item.chargeState(), maxCharge);
      }
    }
  } catch (ifstream::failure e) {
    cerr << "Exception opening/reading file :" << fn << endl;
  } catch (const xercesc::DOMException& e) {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);

  } catch (const xml_schema::exception& e) {
    cerr << e << endl;
  } catch (std::exception e) {
    cerr << e.what() << endl;
  }
  ifs.close();
  return;
}

void MzidentmlReader::read(const std::string &fn, bool isDecoy, boost::shared_ptr<FragSpectrumScanDatabase> database) {
  namespace xml = xsd::cxx::xml;
  scanNumberMapType scanNumberMap;

  try {
    ifstream ifs;
    ifs.exceptions(ifstream::badbit | ifstream::failbit);
    ifs.open(fn.c_str());
    parser p;
    string schemaDefinition = MZIDENTML_SCHEMA_LOCATION + string("mzIdentML1.1.0.xsd");
    string scheme_namespace = MZIDENTML_NAMESPACE;
    string schema_major = "";
    string schema_minor = "";
    xml_schema::dom::auto_ptr<DOMDocument> doc
            (p.start(ifs, fn.c_str(), true, schemaDefinition, schema_major, schema_minor, scheme_namespace));

    //NOTE wouldnt be  better to use the get tag by Name to jump SequenceCollenction directly?
    while (doc.get() != 0 && !XMLString::equals(sequenceCollectionStr,
            doc->getDocumentElement()->getTagName())) {
      doc = p.next(); // Let's skip some sub trees that we are not interested, e.g. AuditCollection
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

    for (doc = p.next(); doc.get() != 0 && !XMLString::equals(spectrumIdentificationResultStr,
            doc->getDocumentElement()->getTagName()); doc = p.next()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }

    unsigned scanNumber = 0;
    for (; doc.get() != 0 && XMLString::equals(spectrumIdentificationResultStr,
            doc->getDocumentElement()->getTagName()); doc = p.next()) {
      ::mzIdentML_ns::SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement());
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge
              = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();

      BOOST_FOREACH(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem()) {
        createPSM(item, experimentalMassToCharge, isDecoy, ++scanNumber, database);
      }

    }

    cleanHashMaps();
    ifs.close();
  } catch (const xercesc::DOMException& e) {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;
    XMLString::release(&tmpStr);
    exit(-1);
  } catch (const xml_schema::exception& e) {
    cerr << e << endl;
    exit(-1);
  } catch (const ios_base::failure&) {
    cerr << "io failure" << endl;
    exit(-1);
  }

  return;
}

void MzidentmlReader::createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
        ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge,
        bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database) {

  std::auto_ptr< percolatorInNs::features > features_p(new percolatorInNs::features());
  percolatorInNs::features::feature_sequence & f_seq = features_p->feature();

  if (!item.calculatedMassToCharge().present()) {
    std::cerr << "error: calculatedMassToCharge attribute is needed for percolator" << std::endl;
    exit(-1);
  }

  std::string peptideSeq = peptideMap[item.peptide_ref().get()]->PeptideSequence();
  std::string peptideId = item.peptide_ref().get();
  std::vector< std::string > proteinIds;
  std::string __flankN = "";
  std::string __flankC = "";

  //FIXME IMPORTANT fix, here I take only 1 peptide per PSM but the option -m might tell me to take more,
  //FIXME I have to modify this loop to obtain more PSMs in that case
  //NOTE I might be able to get the PeptideEVidence and the protein out of the PSM object
  //Get rid of unprintables in proteinName?

  BOOST_FOREACH(const ::mzIdentML_ns::PeptideEvidenceRefType &pepEv_ref, item.PeptideEvidenceRef()) {
    std::string ref_id = pepEv_ref.peptideEvidence_ref().c_str();
    ::mzIdentML_ns::PeptideEvidenceType *pepEv = peptideEvidenceMap[ref_id];
    __flankN = boost::lexical_cast<string > (pepEv->pre());
    __flankC = boost::lexical_cast<string > (pepEv->post());
    if (__flankN == "?") {__flankN = "-";} //MSGF+ sometimes outputs questionmarks here
    if (__flankC == "?") {__flankC = "-";}
    std::string proteinid = boost::lexical_cast<string > (pepEv->dBSequence_ref());
    mzIdentML_ns::SequenceCollectionType::DBSequence_type *proteinObj = proteinMap[proteinid];
    std::string proteinName = boost::lexical_cast<string > (proteinObj->accession());
    proteinIds.push_back(proteinName);
  }

  if (po->iscombined && !po->reversedFeaturePattern.empty()) {
    //NOTE taking the highest ranked PSM protein for combined search
    isDecoy = proteinIds.front().find(po->reversedFeaturePattern, 0) != std::string::npos;
  }
  
  double rank = item.rank();
  double PI = boost::lexical_cast<double>(item.calculatedPI().get());
  int charge = item.chargeState();
  double theoretic_mass = boost::lexical_cast<double>(item.calculatedMassToCharge());
  double observed_mass = boost::lexical_cast<double>(item.experimentalMassToCharge());
  std::string peptideSeqWithFlanks = __flankN + std::string(".") + peptideSeq + std::string(".") + __flankC;
  unsigned peptide_length = peptideLength(peptideSeqWithFlanks);
  double dM = massDiff(observed_mass, theoretic_mass, charge);
  std::map<char, int> ptmMap = po->ptmScheme;
  std::string psmid = boost::lexical_cast<string > (item.id()) + "_" + boost::lexical_cast<string > (useScanNumber) + "_" +
            boost::lexical_cast<string > (charge) + "_" + boost::lexical_cast<string > (rank);
  
  //----------------
  //----if SEQUEST
  //----------------
  if (inputFormat == sequest)
  {
    double lnrSP = 0.0;
    double deltaCN = 0.0;
    double xCorr = 0.0;
    double Sp = 0.0;
    double ionMatched = 0.0;
    double ionTotal = 0.0;

    BOOST_FOREACH(const ::mzIdentML_ns::CVParamType & cv, item.cvParam()) 
    {
      if (cv.value().present()) 
      {
        std::string param_name(cv.name().c_str());
        if (sequestFeatures.count(param_name)) 
        {
          switch (sequestFeatures.at(param_name)) 
          {
            case 0: lnrSP = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 1: deltaCN = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 2: xCorr = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 3: Sp = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 4: ionMatched = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 5: ionTotal = boost::lexical_cast<double>(cv.value().get().c_str());break;
          }
        } 
        else 
        {
          std::cerr << "ERROR : an unmapped Sequest parameter " << param_name << " was not found." << std::endl;
        }
      }
    }
    f_seq.push_back(log(max(1.0, lnrSP)));
    f_seq.push_back(deltaCN);
    f_seq.push_back(xCorr);
    f_seq.push_back(Sp);
    f_seq.push_back(ionMatched / ionTotal);
    f_seq.push_back(observed_mass);
    f_seq.push_back(peptideLength(peptideSeqWithFlanks));
  }
  //----------------
  //----if MS-GF PLUS
  //----------------
  else if (inputFormat == msgfplus)  
  {
    double RawScore = 0.0;
    double DeNovoScore = 0.0;
    double SpecEValue = 0.0;
    double EValue = 0.0;
    int IsotopeError = 0.0;
    double ExplainedIonCurrentRatio = 0.0;
    double NTermIonCurrentRatio = 0.0;
    double CTermIonCurrentRatio = 0.0;
    double MS2IonCurrent = 0.0;

    //Read through cvParam elements
    BOOST_FOREACH(const ::mzIdentML_ns::CVParamType & cv, item.cvParam()) 
    {
      if (cv.value().present()) 
      {
        std::string param_name(cv.name().c_str());
        if (msgfplusFeatures.count(param_name)) 
        {
          switch (msgfplusFeatures.at(param_name)) 
          {
            case 0: RawScore = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 1: DeNovoScore = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 2: SpecEValue = boost::lexical_cast<double>(cv.value().get().c_str());break;
            case 3: EValue = boost::lexical_cast<double>(cv.value().get().c_str());break;
          }
        } 
        else 
        {
          std::cerr << "ERROR : an unmapped MS-GF+ parameter " << param_name << " was not found." << std::endl;
        }
      }
    }
    
    //Read through userParam elements
    BOOST_FOREACH(const ::mzIdentML_ns::UserParamType & up, item.userParam()) 
    {
      if (up.value().present()) 
      {
        std::string param_name(up.name().c_str());
        if (msgfplusFeatures.count(param_name)) 
        {
          switch (msgfplusFeatures.at(param_name)) 
          {
            case 4: IsotopeError = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 5: ExplainedIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str());break;
            case 6: NTermIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str());break;
            case 7: CTermIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str());break;
            case 8: MS2IonCurrent = boost::lexical_cast<double>(up.value().get().c_str());break;
          }
        }
        /*
        else 
        {
          std::cerr << "ERROR : an unmapped MS-GF+ parameter " << param_name << " was not found." << std::endl;
        }*/
      }
    }
    //Add +1 to some features to avoid log(0)
    f_seq.push_back(RawScore);
    f_seq.push_back(log(DeNovoScore+1));
    f_seq.push_back(SpecEValue);
    f_seq.push_back(EValue);
    f_seq.push_back(IsotopeError);
    f_seq.push_back(log(ExplainedIonCurrentRatio+1));
    f_seq.push_back(log(NTermIonCurrentRatio+1));
    f_seq.push_back(log(CTermIonCurrentRatio+1));
    f_seq.push_back(log(MS2IonCurrent));
  }

    f_seq.push_back(observed_mass);
    f_seq.push_back(peptideLength(peptideSeqWithFlanks));
  
    for (int c = minCharge; c <= maxCharge; c++) {
      f_seq.push_back(charge == c ? 1.0 : 0.0); // Charge
    }
    if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) {
      f_seq.push_back(Enzyme::isEnzymatic(peptideSeqWithFlanks.at(0), peptideSeqWithFlanks.at(2)) ? 1.0 : 0.0);
      f_seq.push_back(Enzyme::isEnzymatic(peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 3), peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 1)) ? 1.0 : 0.0);
      f_seq.push_back((double) Enzyme::countEnzymatic(peptideSeq));
    }

    f_seq.push_back(dM);
    f_seq.push_back(abs(dM));

    if (po->calcPTMs) {
      f_seq.push_back(cntPTMs(peptideSeqWithFlanks));
    }
    if (po->pngasef) {
      f_seq.push_back(isPngasef(peptideSeqWithFlanks, isDecoy));
    }
    if (po->calcAAFrequencies) {
      computeAAFrequencies(peptideSeqWithFlanks, f_seq);
    }

    percolatorInNs::occurence::flankN_type flankN = peptideSeqWithFlanks.substr(0, 1);
    percolatorInNs::occurence::flankC_type flankC = peptideSeqWithFlanks.substr(peptideSeqWithFlanks.size() - 1, 1);

    // Strip peptide from termini and modifications 
    std::string peptideS = peptideSeq;
    for (unsigned int ix = 0; ix < peptideSeq.size(); ++ix) {
      if (aaAlphabet.find(peptideSeq[ix]) == string::npos &&
              ambiguousAA.find(peptideSeq[ix]) == string::npos &&
              additionalAA.find(peptideSeq[ix]) == string::npos) {
        if (ptmMap.count(peptideSeq[ix]) == 0) {
          std::cerr << "Peptide sequence " << peptideSeqWithFlanks
                  << " contains modification " << peptideSeq[ix] << " that is not specified by a \"-p\" argument" << std::endl;
          exit(-1);
        }
        peptideSeq.erase(ix, 1);
      }
    }

    std::auto_ptr< percolatorInNs::peptideType > peptide_p(new percolatorInNs::peptideType(peptideSeq));
    // Register the ptms
    for (unsigned int ix = 0; ix < peptideS.size(); ++ix) {
      if (aaAlphabet.find(peptideS[ix]) == string::npos &&
              ambiguousAA.find(peptideS[ix]) == string::npos &&
              additionalAA.find(peptideS[ix]) == string::npos) {
        int accession = ptmMap[peptideS[ix]];
        std::auto_ptr< percolatorInNs::uniMod > um_p(new percolatorInNs::uniMod(accession));
        std::auto_ptr< percolatorInNs::modificationType > mod_p(new percolatorInNs::modificationType(um_p, ix));
        peptide_p->modification().push_back(mod_p);
        peptideS.erase(ix, 1);
      }
    }

    ::percolatorInNs::peptideSpectrumMatch* tmp_psm = new ::percolatorInNs::peptideSpectrumMatch
            (features_p, peptide_p, psmid, isDecoy, observed_mass, theoretic_mass, charge);
    std::auto_ptr< ::percolatorInNs::peptideSpectrumMatch > psm_p(tmp_psm);

    for (std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i) {
      std::auto_ptr< percolatorInNs::occurence > oc_p(new percolatorInNs::occurence(*i, flankN, flankC));
      psm_p->occurence().push_back(oc_p);
    }

    database->savePsm(useScanNumber, psm_p);
    return;
  }
