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

#include "Version.h"

#ifdef XML_SUPPORT

/** some constant strings to be used to compare xml strings **/

// databases
static const XMLCh databasesStr[] = {
    xercesc::chLatin_d, xercesc::chLatin_a, xercesc::chLatin_t, xercesc::chLatin_a,
    xercesc::chLatin_b, xercesc::chLatin_a, xercesc::chLatin_s, xercesc::chLatin_e,
    xercesc::chLatin_s, xercesc::chNull};

// calibration
static const XMLCh calibrationStr[] = {
    xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_l, xercesc::chLatin_i,
    xercesc::chLatin_b, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_t,
    xercesc::chLatin_i, xercesc::chLatin_o, xercesc::chLatin_n, xercesc::chNull};

// proteins
static const XMLCh proteinsStr[] = {
    xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t,
    xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chLatin_s,
    xercesc::chNull};

// protein
static const XMLCh proteinStr[] = {
    xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t,
    xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chNull};

// fragSpectrumScan
static const XMLCh fragSpectrumScanStr[] = {
    xercesc::chLatin_f, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_g,
    xercesc::chLatin_S, xercesc::chLatin_p, xercesc::chLatin_e, xercesc::chLatin_c,
    xercesc::chLatin_t, xercesc::chLatin_r, xercesc::chLatin_u, xercesc::chLatin_m,
    xercesc::chLatin_S, xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_n,
    xercesc::chNull};

#endif  // XML_SUPPORT

XMLInterface::XMLInterface(const std::string& outputFN,
                           const std::string& PEPoutputFN,
                           bool schemaValidation, bool printDecoys, bool printExpMass) : xmlOutputFN_(outputFN), pepXMLOutputFN_(PEPoutputFN), schemaValidation_(schemaValidation), otherCall_(""), reportUniquePeptides_(false), reportPepXML_(false), printDecoys_(printDecoys), printExpMass_(printExpMass) {}

XMLInterface::~XMLInterface() {
    // clean up temporary files if an exception occurred during the writing
    remove(xmlOutputFN_PSMs.c_str());
    remove(xmlOutputFN_Peptides.c_str());
    remove(xmlOutputFN_Proteins.c_str());
}

int XMLInterface::readPin(istream& dataStream, const std::string& xmlInputFN,
                          SetHandler& setHandler, SanityCheck*& pCheck,
                          ProteinProbEstimator* protEstimator, Enzyme*& enzyme) {
    std::vector<double> noWeights;
    Scores noScores(true);
    return readAndScorePin(dataStream, noWeights, noScores, xmlInputFN, setHandler,
                           pCheck, protEstimator, enzyme);
}

int XMLInterface::readAndScorePin(istream& dataStream, std::vector<double>& rawWeights,
                                  Scores& allScores, const std::string& xmlInputFN,
                                  SetHandler& setHandler, SanityCheck*& pCheck,
                                  ProteinProbEstimator* protEstimator, Enzyme*& enzyme) {
#ifdef XML_SUPPORT
    xercesc::XMLPlatformUtils::Initialize();
    try {
        using namespace xercesc;
        string schemaDefinition = Globals::getInstance()->getXMLDir() + PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
        parser p;
        auto deleter = xsd::cxx::xml::dom::deleter<xercesc::DOMDocument>();

        // Create the unique_ptr with the custom deleter
        std::unique_ptr<xercesc::DOMDocument, decltype(deleter)> doc(nullptr, deleter);

        // Initialize the unique_ptr with the document
        doc.reset(p.start(dataStream, xmlInputFN.c_str(), schemaValidation_,
                          schemaDefinition, PIN_VERSION_MAJOR, PIN_VERSION_MINOR).release());

        // Assign new value to doc while maintaining the custom deleter
        doc.reset(p.next().release());

        // read enzyme element
        char* value = XMLString::transcode(doc->getDocumentElement()->getTextContent());

        if (VERB > 1) std::cerr << "enzyme=" << value << std::endl;

        delete enzyme;
        enzyme = Enzyme::createEnzyme(value);

        XMLString::release(&value);
        doc.reset(p.next().release());

        // checking if database is present to jump it
        bool hasProteins = false;
        if (XMLString::equals(databasesStr, doc->getDocumentElement()->getTagName())) {
            doc.reset(p.next().release());
            hasProteins = true;
        }

        // read process_info element
        percolatorInNs::process_info processInfo(*doc->getDocumentElement());
        otherCall_ = processInfo.command_line();
        doc.reset(p.next().release());

        if (XMLString::equals(calibrationStr, doc->getDocumentElement()->getTagName())) {
            doc.reset(p.next().release());
        };

        // read feature names and initial values that are present in feature descriptions
        FeatureNames& featureNames = DataSet::getFeatureNames();
        percolatorInNs::featureDescriptions featureDescriptions(*doc->getDocumentElement());
        for (const auto& featureDesc : featureDescriptions.featureDescription()) {
            featureNames.insertFeature(featureDesc.name());
        }
        featureNames.initFeatures();

        std::vector<double> init_values(FeatureNames::getNumFeatures());
        setHandler.getFeaturePool().createPool(FeatureNames::getNumFeatures());

        bool hasDefaultValues = false;
        unsigned int i = 0;
        for (const auto& featureDesc : featureDescriptions.featureDescription()) {
            if (featureDesc.initialValue().present()) {
                if (featureDesc.initialValue().get() != 0.0)
                    hasDefaultValues = true;
                if (VERB > 2) {
                    std::cerr << "Initial direction for " << featureDesc.name() << " is " << featureDesc.initialValue().get() << std::endl;
                }
                init_values[i] = featureDesc.initialValue().get();
            }
            ++i;
        }

        bool readProteins = true;
        if (rawWeights.empty()) {
            auto targetSet = std::make_unique<DataSet>();
            targetSet->setLabel(LabelType::TARGET);
            auto decoySet = std::make_unique<DataSet>();
            decoySet->setLabel(LabelType::DECOY);

            bool concatenatedSearch = true;

            if (setHandler.getMaxPSMs() > 0u) {
                readProteins = false;
                std::priority_queue<PSMDescriptionPriority> subsetPSMs;
                std::map<ScanId, std::pair<size_t, bool>> scanIdLookUp;
                unsigned int upperLimit = UINT_MAX;
                for (doc.reset(p.next().release());
                     doc && XMLString::equals(fragSpectrumScanStr, doc->getDocumentElement()->getTagName());
                     doc.reset(p.next().release())) {
                    percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
                    for (const auto& psm : fragSpectrumScan.peptideSpectrumMatch()) {
                        ScanId scanId = getScanId(psm, fragSpectrumScan.scanNumber());
                        size_t randIdx;
                        if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
                            if (concatenatedSearch && psm.isDecoy() != scanIdLookUp[scanId].second) {
                                concatenatedSearch = false;
                            }
                            randIdx = scanIdLookUp[scanId].first;
                        } else {
                            randIdx = PseudoRandom::lcg_rand();
                            scanIdLookUp[scanId].first = randIdx;
                            scanIdLookUp[scanId].second = psm.isDecoy();
                        }

                        if (subsetPSMs.size() < setHandler.getMaxPSMs() || randIdx < upperLimit) {
                            PSMDescriptionPriority psmPriority;
                            psmPriority.psm = readPsm(psm, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());
                            psmPriority.label = (psm.isDecoy() ? LabelType::DECOY : LabelType::TARGET);
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

                setHandler.addQueueToSets(subsetPSMs, targetSet.get(), decoySet.get());
            } else {
                std::map<ScanId, bool> scanIdLookUp;
                for (doc.reset(p.next().release());
                     doc && XMLString::equals(fragSpectrumScanStr, doc->getDocumentElement()->getTagName());
                     doc.reset(p.next().release())) {
                    percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
                    for (const auto& psm : fragSpectrumScan.peptideSpectrumMatch()) {
                        ScanId scanId = getScanId(psm, fragSpectrumScan.scanNumber());
                        if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
                            if (concatenatedSearch && psm.isDecoy() != scanIdLookUp[scanId]) {
                                concatenatedSearch = false;
                            }
                        } else {
                            scanIdLookUp[scanId] = psm.isDecoy();
                        }

                        PSMDescription* psmPtr = readPsm(psm, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());
                        if (psm.isDecoy()) {
                            decoySet->registerPsm(psmPtr);
                        } else {
                            targetSet->registerPsm(psmPtr);
                        }
                    }
                }
            }

            setHandler.push_back_dataset(targetSet.release());
            setHandler.push_back_dataset(decoySet.release());

            pCheck = SanityCheck::initialize(otherCall_);
            assert(pCheck);
            pCheck->checkAndSetDefaultDir();
            if (hasDefaultValues) pCheck->addDefaultWeights(init_values);
            pCheck->setConcatenatedSearch(concatenatedSearch);

        } else {
            const unsigned int numFeatures = FeatureNames::getNumFeatures();
            for (doc.reset(p.next().release());
                 doc && XMLString::equals(fragSpectrumScanStr, doc->getDocumentElement()->getTagName());
                 doc.reset(p.next().release())) {
                percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
                for (const auto& psm : fragSpectrumScan.peptideSpectrumMatch()) {
                    ScoreHolder sh;
                    sh.label = (psm.isDecoy() ? LabelType::DECOY : LabelType::TARGET);
                    sh.pPSM = readPsm(psm, fragSpectrumScan.scanNumber(), readProteins, setHandler.getFeaturePool());

                    allScores.scoreAndAddPSM(sh, rawWeights, setHandler.getFeaturePool());
                }
            }
        }

        if (readProteins) {
            unsigned numProteins = 0;
            if (hasProteins && ProteinProbEstimator::getCalcProteinLevelProb()) {
                assert(protEstimator);
                for (doc.reset(p.next().release());
                     doc && XMLString::equals(proteinStr, doc->getDocumentElement()->getTagName());
                     doc.reset(p.next().release())) {
                    ::percolatorInNs::protein protein(*doc->getDocumentElement());
                    protEstimator->addProteinDb(protein.isDecoy(), protein.name(), protein.sequence(), protein.length());
                    ++numProteins;
                }
            }
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
#else   // XML_SUPPORT
    std::cerr << "ERROR: Compiler flag XML_SUPPORT was off, you cannot use the -k flag for pin-format input files" << std::endl;
    return 0;
#endif  // XML_SUPPORT
}

#ifdef XML_SUPPORT
// Convert a peptide with or without modifications into a string
std::string XMLInterface::decoratePeptide(const ::percolatorInNs::peptideType& peptide) {
    std::list<std::pair<int, std::string> > mods;
    std::string peptideSeq = peptide.peptideSequence();
    percolatorInNs::peptideType::modification_const_iterator modIt;
    modIt = peptide.modification().begin();
    for (; modIt != peptide.modification().end(); ++modIt) {
        std::stringstream ss;
        if (modIt->uniMod().present()) {
            ss << "[UNIMOD:" << modIt->uniMod().get().accession() << "]";
            mods.push_back(std::pair<int, std::string>(modIt->location(), ss.str()));
        }
        if (modIt->freeMod().present()) {
            ss << "[" << modIt->freeMod().get().moniker() << "]";
            mods.push_back(std::pair<int, std::string>(modIt->location(), ss.str()));
        }
    }
    mods.sort(greater<std::pair<int, std::string> >());
    std::list<std::pair<int, std::string> >::const_iterator it;
    for (it = mods.begin(); it != mods.end(); ++it) {
        peptideSeq.insert(it->first, it->second);
    }
    return peptideSeq;
}

PSMDescription* XMLInterface::readPsm(
    const percolatorInNs::peptideSpectrumMatch& psm, unsigned scanNumber,
    bool readProteins, FeatureMemoryPool& featurePool) {
    PSMDescription* myPsm = new PSMDescription();
    string mypept = decoratePeptide(psm.peptide());

    if (psm.occurence().size() <= 0) {
        ostringstream temp;
        temp << "Error: cannot add PSM " << psm.id() << " to the dataset.\n\
    The PSM does not contain protein occurences."
             << std::endl;
        throw MyException(temp.str());
    }

    percolatorInNs::peptideSpectrumMatch::occurence_const_iterator occIt;
    occIt = psm.occurence().begin();
    for (; occIt != psm.occurence().end(); ++occIt) {
        if (readProteins) myPsm->proteinIds.push_back(occIt->proteinId());
        // adding n-term and c-term residues to peptide
        // NOTE the residues for the peptide in the PSMs are always the same for every protein
        myPsm->setPeptide(occIt->flankN() + "." + mypept + "." + occIt->flankC());
    }

    myPsm->setId(psm.id());
    myPsm->scan = scanNumber;
    myPsm->expMass = psm.experimentalMass();
    myPsm->calcMass = psm.calculatedMass();
    if (psm.observedTime().present()) {
        myPsm->setRetentionTime(psm.observedTime().get());
    }

    myPsm->features = featurePool.allocate();

    for (unsigned int i = 0; i < psm.features().feature().size(); ++i) {
        myPsm->features[i] = psm.features().feature()[i];
    }

    return myPsm;
}

ScanId XMLInterface::getScanId(const percolatorInNs::peptideSpectrumMatch& psm,
                               unsigned scanNumber) {
    ScanId scanId;
    scanId.first = scanNumber;
    scanId.second = psm.experimentalMass();
    return scanId;
}
#endif  // XML_SUPPORT

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
    os << "  </psms>" << endl
       << endl;
    os.close();
}

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
    os << "  </peptides>" << endl
       << endl;
    os.close();
}

/**
 * Subroutine of @see XMLInterface::writeXML() for protein output
 */
void XMLInterface::writeXML_Proteins(ProteinProbEstimator* protEstimator) {
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
       << endl
       << "xmlns=\"" << space << "\" "
       << endl
       << "xmlns:p=\"" << space << "\" "
       << endl
       << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
       << endl
       << "xsi:schemaLocation=\"" << schema << "\" "
       << endl
       << "p:majorVersion=\"" << VERSION_MAJOR << "\" p:minorVersion=\""
       << VERSION_MINOR << "\" p:percolator_version=\"Percolator version "
       << VERSION << "\">\n"
       << endl;
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
    os << "    <psms_qlevel>" << numberQpsms_ << "</psms_qlevel>" << endl;
    if (reportUniquePeptides_)
        os << "    <peptides_qlevel>" << fullset.getQvaluesBelowLevel(0.01) << "</peptides_qlevel>" << endl;
    if (ProteinProbEstimator::getCalcProteinLevelProb())
        os << "    <proteins_qlevel>" << protEstimator->getQvaluesBelowLevel(0.01) << "</proteins_qlevel>" << endl;
    os << "  </process_info>" << endl
       << endl;

    // append PSMs
    ifstream ifs_psms(xmlOutputFN_PSMs.data(), ios::in | ios::binary);
    os << ifs_psms.rdbuf();
    ifs_psms.close();
    remove(xmlOutputFN_PSMs.c_str());
    // append Peptides
    if (reportUniquePeptides_) {
        ifstream ifs_peptides(xmlOutputFN_Peptides.data(), ios::in | ios::binary);
        os << ifs_peptides.rdbuf();
        ifs_peptides.close();
        remove(xmlOutputFN_Peptides.c_str());
    }
    // append Proteins
    if (ProteinProbEstimator::getCalcProteinLevelProb()) {
        ifstream ifs_proteins(xmlOutputFN_Proteins.data(), ios::in | ios::binary);
        os << ifs_proteins.rdbuf();
        ifs_proteins.close();
        remove(xmlOutputFN_Proteins.c_str());
    }
    os << "</percolator_output>" << endl;
    os.close();
}

void XMLInterface::writePepXML(Scores& fullset, ProteinProbEstimator* protEstimator, std::string call) {
    ofstream os;
    const string schema =  // space +
        "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v122.xsd";

    os.open(pepXMLOutputFN_.data(), ios::out | ios::binary);
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    os << "<msms_pipeline_analysis date=\"" << getAtomicTime() << "\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" summary_xml=\"/sdd/proteomics/DATA/PXD014076/iProphet3/interact-Symb_Proteome_DIA_RAW_S05_Q3.pep.xml\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl;
    os << "    <analysis_summary analysis=\"peptideprophet\" time=\"" << getAtomicTime() << "\">" << endl;
    os << "    </analysis_summary>" << endl;

    std::string pepPath = pepXMLOutputFN_.data();

    /* TODO: MUST FIX */
    std::string base_name = pepPath.substr(0, pepPath.find("."));

    ifstream ifs_psms(xmlpeptideOutputFN_PSMs.data(), ios::in | ios::binary);
    os << ifs_psms.rdbuf();
    ifs_psms.close();
    remove(xmlpeptideOutputFN_PSMs.c_str());

    os << "    </msms_run_summary>" << endl;

    os << "</msms_pipeline_analysis>" << endl;
    os.close();
}

/* Get time-stamp for pepXML */
void XMLInterface::setAtomicTime() {
    time_t now = time(0);

    // convert now to string form
    tm* ltm = localtime(&now);

    std::string date = "";
    std::string time = "";

    int year = ltm->tm_year;
    date += to_string(1900 + year);

    int month = 1 + ltm->tm_mon;
    std::string month_s = to_string(month);

    if (month_s.length() < 2) {
        date += "-0" + month_s;
    } else {
        date += "-" + month_s;
    }

    int day = ltm->tm_mday;
    std::string day_s = to_string(day);

    if (day_s.length() < 2) {
        date += "-0" + day_s;
    } else {
        date += "-" + day_s;
    }

    int hour = 5 + ltm->tm_hour;
    std::string hour_s = to_string(hour);
    if (hour_s.length() < 2) {
        time += "0" + hour_s;
    } else {
        time += hour_s;
    }

    int min = ltm->tm_min;

    std::string min_s = to_string(min);
    if (min_s.length() < 2) {
        time += ":0" + min_s;
    } else {
        time += ":" + min_s;
    }

    int sec = ltm->tm_sec;
    std::string sec_s = to_string(sec);
    if (sec_s.length() < 2) {
        time += ":0" + sec_s;
    } else {
        time += ":" + sec_s;
    }

    atomicDate = date + "T" + time;
}

// Change this function name to something including a string pepXML
void XMLInterface::writePepXML_PSMs(Scores& fullset, double selectionFdr_, std::string protEstimatorDecoyPrefix) {
    setAtomicTime();

    pi0Psms_ = fullset.getPi0();
    numberQpsms_ = fullset.getQvaluesBelowLevel(0.01);

    ofstream os;

    /* xmlpeptideOutputFN_ = pepXMLOutputFN_; */
    xmlpeptideOutputFN_PSMs.append("writePepXML_PSMs");
    os.open(xmlpeptideOutputFN_PSMs.c_str(), ios::out);

    /* os << "  <psms>" << endl; */
    map<char, float> aaDict = getRoughAminoWeightDict();
    /* Sort psms based on base name  */
    std::sort(fullset.begin(), fullset.end(), lessThanBaseName());
    std::string pepXMLBaseName = "";
    bool first_msms_summary = true;

    int index = 1;
    for (std::vector<ScoreHolder>::iterator sh = fullset.begin(); sh != fullset.end(); ++sh) {
        std::string id = sh->pPSM->getId();
        auto baseName = id.substr(0, id.find('.'));
        if (baseName != pepXMLBaseName) {
            /*  Start of a new msms run*/
            pepXMLBaseName = baseName;
            if (first_msms_summary) {
                first_msms_summary = false;
            } else {
                /* End of msms run */
                os << "    </msms_run_summary>" << endl;
            }

            /* New msms run! */
            os << "    <msms_run_summary base_name=\"" << baseName << "\" raw_data_type=\"mzML\" raw_data=\"mzML\">" << endl;
            os << "        <search_summary base_name=\"" << baseName << "\" search_engine=\"X! Tandem\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">" << endl;
            os << "            <parameter name=\"decoy_prefix\" value=\"" << protEstimatorDecoyPrefix << "\" />" << endl;
            os << "        </search_summary>" << endl;
            os << "        <analysis_timestamp analysis=\"peptideprophet\" time=\"" << getAtomicTime() << "\" id=\"1\"/>" << endl;
        }
        if (sh->q < selectionFdr_)
            sh->printPepXML(os, aaDict, index);
        index++;
    }
    /* os << "  </psms>" << endl << endl; */
    os.close();
}

map<char, float> XMLInterface::getRoughAminoWeightDict() {
    // This function gives an approximate amino acid weight
    // with two digits precision. For parsing purposes.
    map<char, float> roughAminoAcidWeight =
        boost::assign::map_list_of('A', 71.04)('C', 103.01)('D', 115.03)('E', 129.04)('F', 147.07)('G', 57.02)('H', 137.06)('I', 113.08)('K', 128.09)('L', 113.08)('M', 131.04)('N', 114.04)('P', 97.05)('Q', 128.06)('R', 156.10)('S', 87.03)('T', 101.05)('V', 99.07)('W', 186.08)('Y', 163.06);
    return roughAminoAcidWeight;
}
