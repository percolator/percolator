#if defined(_MSC_VER)
#define _CMATH_IN_CRT 1
#include <float.h>
#endif

#include <xsd/cxx/xml/dom/auto-ptr.hxx> // For xsd::cxx::xml::dom::deleter

#include "SequestReader.h"


const std::map<string, int> SequestReader::sequestFeatures =
boost::assign::map_list_of("sequest:PeptideRankSp", 0)
                                  ("sequest:deltacn", 1)
                                  ("sequest:xcorr", 2)
                                  ("sequest:PeptideSp", 3)
                                  ("sequest:matched ions", 4)
                                  ("sequest:total ions", 5)
                                  ("sequest:PeptideIdnumber", 6)
                                  ("sequest:PeptideNumber", 7);

SequestReader::SequestReader(ParseOptions po) : MzidentmlReader(po) {}

SequestReader::~SequestReader() {}

bool SequestReader::checkValidity(const std::string &file) {

  bool isvalid = true;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  if (!fileIn) {
    ostringstream temp;
    temp << "Error : can not open file " << file << std::endl;
    isvalid = false;
    throw MyException(temp.str());
  }
  std::string line;
  if (!getline(fileIn, line)) {
    ostringstream temp;
    temp << "Error : can not read file " << file << std::endl;
    isvalid = false;
    throw MyException(temp.str());
  }
  if (line.find("<?xml") == std::string::npos) {
    fileIn.close();
    ostringstream temp;
    temp << "Error : the input file is not xml format " << file << std::endl;
    isvalid = false;
    throw MyException(temp.str());
  }
  else //Test whether Sequest or MS-GF+ format
  {
    std::string line2, line3;
    getline(fileIn, line2);
    getline(fileIn, line3);

    if ((line2[1] != '!' && line2.find("SEQUEST") != std::string::npos && line2.find("MzIdentML") != std::string::npos)
         || (line3[1] != '!' && line3.find("SEQUEST") != std::string::npos && line3.find("MzIdentML") != std::string::npos))
    {
      if(VERB > 2)
	std::cerr << "MzIdentML - SEQUEST format" << std::endl;
      isvalid = true;
    } else {
      fileIn.close();
      ostringstream temp;
      temp << "Error : the input file is not MzIdentML - Sequest format " << file << std::endl;
      isvalid = false;
      throw MyException(temp.str());
    }

  }
  fileIn.close();
  return isvalid;
}



void SequestReader::addFeatureDescriptions(bool doEnzyme)
{

  push_backFeatureDescription("lnrSp");
  push_backFeatureDescription("deltCn");
  push_backFeatureDescription("Xcorr");
  push_backFeatureDescription("Sp");
  push_backFeatureDescription("IonFrac");
  push_backFeatureDescription("Mass");
  push_backFeatureDescription("PepLen");
  push_backFeatureDescription("dM");
  push_backFeatureDescription("absdM");

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

  if (po.calcPTMs) {
    push_backFeatureDescription("ptm");
  }
  if (po.pngasef) {
    push_backFeatureDescription("PNGaseF");
  }

  if (po.calcAAFrequencies) {
    BOOST_FOREACH (const char aa, freqAA) {
      std::string temp = std::string(1,aa) + "-Freq";
      push_backFeatureDescription(temp.c_str());
    }
  }

}

void SequestReader::createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
        ::percolatorInNs::fragSpectrumScan::experimentalMass_type experimentalMass,
        bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database,
        const std::string & fn) {

  std::unique_ptr<percolatorInNs::features> features_p(new percolatorInNs::features());
  percolatorInNs::features::feature_sequence & f_seq = features_p->feature();

  if (!item.calculatedMassToCharge().present()) {
    ostringstream temp;
    temp << "Error: calculatedMassToCharge attribute not found in PSM "
    << boost::lexical_cast<string>(item.id()) << std::endl;
    throw MyException(temp.str());
  }

  std::string peptideSeq = peptideMap[item.peptide_ref().get()]->PeptideSequence();
  std::string peptideId = item.peptide_ref().get();
  std::vector<std::string> proteinIds;
  std::string __flankN = "";
  std::string __flankC = "";

  try {
    BOOST_FOREACH (const ::mzIdentML_ns::PeptideEvidenceRefType &pepEv_ref, item.PeptideEvidenceRef()) {
      std::string ref_id = pepEv_ref.peptideEvidence_ref().c_str();
      ::mzIdentML_ns::PeptideEvidenceType *pepEv = peptideEvidenceMap[ref_id];
      if (peptideId != std::string(pepEv->peptide_ref())) {
        std::cerr << "Warning: The PSM " << boost::lexical_cast<string>(item.id())
                  << " contains different chimeric peptide sequences. "
                  << peptideMap[pepEv->peptide_ref()]->PeptideSequence() << " and " << peptideSeq
                  << " only the proteins that contain the first peptide will be included in the PSM.\n" << std::endl;
      } else {
        __flankN = boost::lexical_cast<string>(pepEv->pre());
        __flankC = boost::lexical_cast<string>(pepEv->post());
        if (__flankN == "?") { __flankN = "-"; }
        if (__flankC == "?") { __flankC = "-"; }
        std::string proteinid = boost::lexical_cast<string>(pepEv->dBSequence_ref());
        mzIdentML_ns::SequenceCollectionType::DBSequence_type *proteinObj = proteinMap[proteinid];
        std::string proteinName = boost::lexical_cast<string>(proteinObj->accession());
        proteinIds.push_back(proteinName);
      }
    }

    if (__flankC.empty() || __flankN.empty()) {
      ostringstream temp;
      temp << "Error: The PSM " << boost::lexical_cast<string>(item.id()) << " is ill-formed." << std::endl;
      throw MyException(temp.str());
    }

    if (po.iscombined && !po.reversedFeaturePattern.empty()) {
      isDecoy = proteinIds.front().find(po.reversedFeaturePattern, 0) != std::string::npos;
    }

    double rank = item.rank();
    double PI = boost::lexical_cast<double>(item.calculatedPI().get());
    int charge = item.chargeState();
    double theoretic_mass = boost::lexical_cast<double>(item.calculatedMassToCharge());
    double observed_mass = boost::lexical_cast<double>(item.experimentalMassToCharge());
    std::string peptideSeqWithFlanks = __flankN + std::string(".") + peptideSeq + std::string(".") + __flankC;
    unsigned peptide_length = peptideLength(peptideSeqWithFlanks);
    std::map<char, int> ptmMap = po.ptmScheme;
    std::string peptideNoMods = removePTMs(peptideSeqWithFlanks, ptmMap);
    std::string psmId = createPsmId(item.id(), observed_mass, useScanNumber, charge, static_cast<unsigned int>(rank));

    double lnrSP = 0.0;
    double deltaCN = 0.0;
    double xCorr = 0.0;
    double Sp = 0.0;
    double ionMatched = 0.0;
    double ionTotal = 0.0;
    double dM = massDiff(observed_mass, theoretic_mass, static_cast<unsigned int>(charge));

    BOOST_FOREACH (const ::mzIdentML_ns::CVParamType & cv, item.cvParam()) {
      if (cv.value().present()) {
        std::string param_name(cv.name().c_str());
        if (sequestFeatures.count(param_name)) {
          switch (sequestFeatures.at(param_name)) {
            case 0: lnrSP = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 1: deltaCN = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 2: xCorr = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 3: Sp = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 4: ionMatched = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 5: ionTotal = boost::lexical_cast<double>(cv.value().get().c_str()); break;
          }
        } else {
          std::cerr << "Error: An unmapped Sequest parameter " << param_name << " was not found." << std::endl;
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
    f_seq.push_back(dM);
    f_seq.push_back(abs(dM));

    for (int c = minCharge; c <= maxCharge; c++) {
      f_seq.push_back(charge == c ? 1.0 : 0.0); // Charge
    }
    if (enzyme_->getEnzymeType() != Enzyme::NO_ENZYME) {
      f_seq.push_back(enzyme_->isEnzymatic(peptideNoMods.at(0), peptideNoMods.at(2)) ? 1.0 : 0.0);
      f_seq.push_back(enzyme_->isEnzymatic(peptideNoMods.at(peptideNoMods.size() - 3), peptideNoMods.at(peptideNoMods.size() - 1)) ? 1.0 : 0.0);
      std::string peptide2 = peptideNoMods.substr(2, peptideNoMods.size() - 4);
      f_seq.push_back(static_cast<double>(enzyme_->countEnzymatic(peptide2)));
    }

    if (po.calcPTMs) {
      f_seq.push_back(cntPTMs(peptideSeqWithFlanks, ptmMap));
    }
    if (po.pngasef) {
      f_seq.push_back(isPngasef(peptideSeqWithFlanks, isDecoy));
    }
    if (po.calcAAFrequencies) {
      computeAAFrequencies(peptideSeqWithFlanks, f_seq);
    }

    percolatorInNs::occurence::flankN_type flankN = peptideSeqWithFlanks.substr(0, 1);
    percolatorInNs::occurence::flankC_type flankC = peptideSeqWithFlanks.substr(peptideSeqWithFlanks.size() - 1, 1);

    std::string peptideS = peptideSeq;
    for (unsigned int ix = 0; ix < peptideSeq.size(); ++ix) {
      if (freqAA.find(peptideSeq[ix]) == string::npos) {
        if (ptmMap.count(peptideSeq[ix]) == 0) {
          ostringstream temp;
          temp << "Error: Peptide sequence " << peptideSeqWithFlanks
               << " contains modification " << peptideSeq[ix] << " that is not specified by a \"-p\" argument" << std::endl;
          throw MyException(temp.str());
        }
        peptideSeq.erase(ix--, 1);
      }
    }

    std::unique_ptr<percolatorInNs::peptideType> peptide_p(new percolatorInNs::peptideType(peptideSeq));
    for (unsigned int ix = 0; ix < peptideS.size(); ++ix) {
      if (freqAA.find(peptideS[ix]) == string::npos) {
        int accession = ptmMap[peptideS[ix]];
        std::unique_ptr<percolatorInNs::uniMod> um_p(new percolatorInNs::uniMod(accession));
        std::unique_ptr<percolatorInNs::modificationType> mod_p(new percolatorInNs::modificationType(static_cast<int>(ix)));
        mod_p->uniMod(std::move(um_p)); // Use std::move to transfer ownership
        peptide_p->modification().push_back(std::move(mod_p)); // Use std::move to transfer ownership
        peptideS.erase(ix--, 1);
      }
    }

    std::unique_ptr<percolatorInNs::peptideSpectrumMatch> psm_p(new percolatorInNs::peptideSpectrumMatch(
      std::move(features_p), std::move(peptide_p), psmId, isDecoy, observed_mass, theoretic_mass, charge));

    for (const std::string &proteinId : proteinIds) {
      std::unique_ptr<percolatorInNs::occurence> oc_p(new percolatorInNs::occurence(proteinId, flankN, flankC));
      psm_p->occurence().push_back(std::move(oc_p)); // Use std::move to transfer ownership
    }

    database->savePsm(useScanNumber, std::move(psm_p)); // Use std::move to transfer ownership
  }
  catch(std::exception const& e) {
    ostringstream temp;
    temp << "Error: parsing PSM: " << boost::lexical_cast<string>(item.id())
         << "\nThe error was: " << e.what() << std::endl;
    throw MyException(temp.str());
  }

  return;
}

