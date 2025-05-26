#if defined(_MSC_VER)
#define _CMATH_IN_CRT 1
#include <float.h>
#endif

#include "MsgfplusReader.h"

//NOTE ugly hack to get the order of the values of the features according to their names
const std::map<string, int> MsgfplusReader::msgfplusFeatures =
  boost::assign::map_list_of("MS-GF:RawScore", 0)
    ("MS-GF:DeNovoScore", 1)
    ("MS-GF:SpecEValue", 2)
    ("MS-GF:EValue", 3)
    //All below features are on user specified element userParam
    ("IsotopeError", 4)
    ("ExplainedIonCurrentRatio", 5)
    ("NTermIonCurrentRatio", 6)
    ("CTermIonCurrentRatio", 7)
    ("MS2IonCurrent", 8)
    ("MeanRelErrorTop7", 9)
    ("StdevRelErrorTop7", 10)
    ("NumMatchedMainIons", 11);

//default score vector //TODO move this to a file or input parameter                          
const std::map<string,double> MsgfplusReader::msgfplusFeaturesDefaultValue =
  boost::assign::map_list_of("RawScore", 0.0)
    ("DeNovoScore",-1.0)
    ("ScoreRatio", 0.0)
    ("Energy", -2.0)
    ("lnEValue", 6.0)
    //("lnSpecEValue", 10.0)
    ("IsotopeError", -2.5)
    ("lnExplainedIonCurrentRatio", 0.0)
    ("lnNTermIonCurrentRatio", 0.0)
    ("lnCTermIonCurrentRatio", 0.0)
    ("lnMS2IonCurrent", 0.0)
    ("Mass", 0.0)
    ("PepLen", 0.0)
    ("dM", 0.0)
    ("absdM", -1.0);

MsgfplusReader::MsgfplusReader(ParseOptions po) :
		MzidentmlReader(po),
		useFragmentSpectrumFeatures(false),
		additionalMsgfFeatures(false),
		numMatchedIonLimit(7),
		neutron(1.0033548378) {
}

MsgfplusReader::~MsgfplusReader() {
}

bool MsgfplusReader::checkValidity(const std::string &file) {

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

    if ((line2[1] != '!' && line2.find("MS-GF+") != std::string::npos && line2.find("MzIdentML") != std::string::npos)
         || (line3[1] != '!' && line3.find("MS-GF+") != std::string::npos && line3.find("MzIdentML") != std::string::npos))
    {
      if(VERB > 2)
        std::cerr << "MzIdentML - MSGF+ format" << std::endl;
      isvalid = true;
    } else {
      fileIn.close();
      ostringstream temp;
      temp << "Error : the input file is not MzIdentML - MSGF+ format " << file << std::endl;
      isvalid = false;
      throw MyException(temp.str());
    }

  }
  fileIn.close();
  return isvalid;
}

void MsgfplusReader::searchEngineSpecificParsing(
	const ::mzIdentML_ns::SpectrumIdentificationItemType & item, const int itemCount) {
	// First, check whether addFeatures was set to 1, in MS-GF+
	if (!additionalMsgfFeatures) {
    	BOOST_FOREACH (const ::mzIdentML_ns::UserParamType & up, item.userParam()) {
    		if (up.value().present()) {
    			std::string param_name(up.name().c_str());
    			// Check whether the mzid-file seem to include the additional features
    			if (param_name == "ExplainedIonCurrentRatio") {  // If one additional feature is found
    				additionalMsgfFeatures = true;
    			}
    		}
    	}
    	if (!additionalMsgfFeatures) {  // If no additional features were found in first PSM
    		ostringstream temp;
    		temp << "Error: No features for learning were found in the mzid-file."
    		<< " Run MS-GF+ with the addFeatures option set to 1." << std::endl;
    		throw MyException(temp.str());
    	}
    }

	// Check whether fragmentation spectrum features are present
	if (!useFragmentSpectrumFeatures) {
  	BOOST_FOREACH(const ::mzIdentML_ns::UserParamType & up, item.userParam()) {
  		if (up.value().present()) {
  			std::string param_name(up.name().c_str());
  			// Check whether the mzid-file seem to include features for fragment spectra resolution and accuracy
  			if (param_name == "MeanRelErrorTop7") {  // If one fragmentSpectrum feature is found
  				useFragmentSpectrumFeatures = true;
  				std::cerr << "Uses features for fragment spectra mass errors" << std::endl;
  			}
  		}
  	}
  }
}


void MsgfplusReader::addFeatureDescriptions(bool doEnzyme)
{

  push_backFeatureDescription("RawScore","",msgfplusFeaturesDefaultValue.at("RawScore"));
  push_backFeatureDescription("DeNovoScore","",msgfplusFeaturesDefaultValue.at("DeNovoScore"));
  push_backFeatureDescription("ScoreRatio","",msgfplusFeaturesDefaultValue.at("ScoreRatio"));
  push_backFeatureDescription("Energy","",msgfplusFeaturesDefaultValue.at("Energy"));
  push_backFeatureDescription("lnEValue","",msgfplusFeaturesDefaultValue.at("lnEValue"));
  //push_backFeatureDescription("lnSpecEValue","",msgfplusFeaturesDefaultValue.at("lnSpecEValue")); // causes problems when used together with lnEValue
  //The below are from element userParam
  push_backFeatureDescription("IsotopeError","",msgfplusFeaturesDefaultValue.at("IsotopeError"));
  push_backFeatureDescription("lnExplainedIonCurrentRatio","",msgfplusFeaturesDefaultValue.at("lnExplainedIonCurrentRatio"));
  push_backFeatureDescription("lnNTermIonCurrentRatio","",msgfplusFeaturesDefaultValue.at("lnNTermIonCurrentRatio"));
  push_backFeatureDescription("lnCTermIonCurrentRatio","",msgfplusFeaturesDefaultValue.at("lnCTermIonCurrentRatio"));
  push_backFeatureDescription("lnMS2IonCurrent","",msgfplusFeaturesDefaultValue.at("lnMS2IonCurrent"));
  push_backFeatureDescription("Mass","",msgfplusFeaturesDefaultValue.at("Mass"));
  push_backFeatureDescription("PepLen","",msgfplusFeaturesDefaultValue.at("PepLen"));
  push_backFeatureDescription("dM","",msgfplusFeaturesDefaultValue.at("dM"));
  push_backFeatureDescription("absdM","",msgfplusFeaturesDefaultValue.at("absdM"));

  //the rest of the features will get default value 0.0
  
  if (useFragmentSpectrumFeatures) {
	  push_backFeatureDescription("MeanErrorTop7");
	  push_backFeatureDescription("sqMeanErrorTop7");
	  push_backFeatureDescription("StdevErrorTop7");
  }

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


double MsgfplusReader::rescaleFragmentFeature(double featureValue, int NumMatchedMainIons) {
	// Rescale the fragment features to penalize features calculated by few ions
	int numerator = (1+numMatchedIonLimit)*(1+numMatchedIonLimit);
	int denominator = (1+(std::min)(NumMatchedMainIons, numMatchedIonLimit))*(1+(std::min)(NumMatchedMainIons, numMatchedIonLimit));
	return featureValue * ((double)numerator/denominator);
}


void MsgfplusReader::createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
        ::percolatorInNs::fragSpectrumScan::experimentalMass_type experimentalMass,
        bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database,
        const std::string &fn) {

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
      // NOTE check that there are not chimeric peptides
      if (peptideId != std::string(pepEv->peptide_ref())) {
        std::cerr << "Warning: The PSM " << boost::lexical_cast<string>(item.id())
                  << " contains different chimeric peptide sequences. "
                  << peptideMap[pepEv->peptide_ref()]->PeptideSequence() << " and " << peptideSeq
                  << " only the proteins that contain the first peptide will be included in the PSM.\n" << std::endl;
      }
      if (__flankN != "-") {
        __flankN = boost::lexical_cast<std::string>(pepEv->pre());
        if (__flankN == "?") { __flankN = "-"; } // MSGF+ sometimes outputs question marks here
        // MT: MSGF+ clips methionine of protein N-terminals, set to "-" to avoid confusion with cleavage rules
        if (__flankN == "M" && boost::lexical_cast<std::string>(pepEv->start()) == "2") { __flankN = "-"; }
      }
      
      if (__flankC != "-") {
        __flankC = boost::lexical_cast<std::string>(pepEv->post());
        if (__flankC == "?") { __flankC = "-"; }
      }
      
      std::string proteinid = boost::lexical_cast<string>(pepEv->dBSequence_ref());
      mzIdentML_ns::SequenceCollectionType::DBSequence_type *proteinObj = proteinMap[proteinid];
      std::string proteinName = boost::lexical_cast<string>(proteinObj->accession());
      proteinIds.push_back(proteinName);
    }

    if (__flankC.empty() || __flankN.empty()) {
      ostringstream temp;
      temp << "Error: The PSM " << boost::lexical_cast<string>(item.id()) << " is ill-formed." << std::endl;
      throw MyException(temp.str());
    }

    if (po.iscombined && !po.reversedFeaturePattern.empty()) {
      // NOTE taking the highest ranked PSM protein for combined search
      isDecoy = proteinIds.front().find(po.reversedFeaturePattern, 0) != std::string::npos;
    }

    double rank = item.rank();
    int charge = item.chargeState();
    double theoretic_mass = boost::lexical_cast<double>(item.calculatedMassToCharge().get());
    double observed_mass = boost::lexical_cast<double>(item.experimentalMassToCharge());
    
    std::string peptideSeqWithFlanks = __flankN + std::string(".") + peptideSeq + std::string(".") + __flankC;
    unsigned peptide_length = peptideLength(peptideSeqWithFlanks);

    // Make a PSM id, from filename, item_id, scan_number, charge and rank
    std::string fileId = fn;
    size_t spos = fileId.rfind('/');
    if (spos != std::string::npos) {
      fileId.erase(0, spos + 1);
    }
    spos = fileId.find('.');
    if (spos != std::string::npos) {
      fileId.erase(spos);
    }
    std::string psmId = createPsmId(fileId + "_" + boost::lexical_cast<string>(item.id()), 
        observed_mass, useScanNumber, charge, static_cast<unsigned int>(rank));

    double RawScore = 0.0;
    double DeNovoScore = 0.0;
    double SpecEValue = 0.0;
    double EValue = 0.0;
    int IsotopeError = 0;
    double ExplainedIonCurrentRatio = 0.0;
    double NTermIonCurrentRatio = 0.0;
    double CTermIonCurrentRatio = 0.0;
    double MS2IonCurrent = 0.0;
    // fragmentFeatureValues
    double MeanErrorTop7 = 0.0;
    double StdevErrorTop7 = 0.0;
    int NumMatchedMainIons = 0;

    // Read through cvParam elements
    BOOST_FOREACH(const ::mzIdentML_ns::CVParamType & cv, item.cvParam()) {
      if (cv.value().present()) {
        std::string param_name(cv.name().c_str());
        if (msgfplusFeatures.count(param_name)) {
          switch (msgfplusFeatures.at(param_name)) {
            case 0: RawScore = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 1: DeNovoScore = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 2: SpecEValue = boost::lexical_cast<double>(cv.value().get().c_str()); break;
            case 3: EValue = boost::lexical_cast<double>(cv.value().get().c_str()); break;
          }
        }
      }
    }

    // Read through userParam elements
    BOOST_FOREACH(const ::mzIdentML_ns::UserParamType & up, item.userParam()) {
      // If a feature has a value NaN, the default values from initialization are used
      if (up.value().present() && boost::lexical_cast<string>(up.value().get().c_str()) != "NaN") {
        std::string param_name(up.name().c_str());
        if (msgfplusFeatures.count(param_name)) {
          switch (msgfplusFeatures.at(param_name)) {
            case 4: IsotopeError = boost::lexical_cast<int>(up.value().get().c_str()); break;
            case 5: ExplainedIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 6: NTermIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 7: CTermIonCurrentRatio = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 8: MS2IonCurrent = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 9: MeanErrorTop7 = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 10:
              // Stdev could equal 0, use the mean error in that case
              if (boost::lexical_cast<string>(up.value().get().c_str()) == "0.0") StdevErrorTop7 = MeanErrorTop7;
              else StdevErrorTop7 = boost::lexical_cast<double>(up.value().get().c_str()); break;
            case 11: NumMatchedMainIons = boost::lexical_cast<int>(up.value().get().c_str()); break;
          }
        }
      } else {
        std::cerr << "PSM: " << boost::lexical_cast<string>(item.id()) << " has feature with value NaN, "
                  << "use the default value for that feature." << std::endl;
      }
    }

    // If MeanErrorAll is 0.0, it was not updated, it was probably missing in the file.
    if (useFragmentSpectrumFeatures && MeanErrorTop7 == 0.0) {
      // Skip this PSM
      // std::cerr << "Skipping PSM with id " << psmId << " because MeanErrorTop7 = 0" << std::endl; // disabled this warning because it occurs a lot
      return;
    }
    
    // The raw theoretical mass from MSGF+ is often of the wrong isotope
    double dM = (observed_mass - (IsotopeError * neutron / charge) - theoretic_mass) / observed_mass;
    // double dM = massDiff(observed_mass, theoretic_mass, charge);  // Gives trouble because of isotopes
    // Add a small number to some logged features to avoid log(0)
    f_seq.push_back(RawScore);
    f_seq.push_back(DeNovoScore);
    f_seq.push_back((std::max)(-1.0, RawScore / (DeNovoScore + 0.0001)));  // ScoreRatio
    f_seq.push_back(DeNovoScore - RawScore);  // Score difference, or Energy
    f_seq.push_back(-log(EValue));
    f_seq.push_back(IsotopeError);
    f_seq.push_back(log(ExplainedIonCurrentRatio + 0.0001));
    f_seq.push_back(log(NTermIonCurrentRatio + 0.0001));
    f_seq.push_back(log(CTermIonCurrentRatio + 0.0001));
    f_seq.push_back(log(MS2IonCurrent));
    f_seq.push_back(observed_mass);
    f_seq.push_back(peptideLength(peptideSeqWithFlanks));
    f_seq.push_back(dM);
    f_seq.push_back(abs(dM));

    if (useFragmentSpectrumFeatures) {
      f_seq.push_back(rescaleFragmentFeature(MeanErrorTop7, NumMatchedMainIons));
      f_seq.push_back(rescaleFragmentFeature(MeanErrorTop7 * MeanErrorTop7, NumMatchedMainIons));  // squared
      f_seq.push_back(rescaleFragmentFeature(StdevErrorTop7, NumMatchedMainIons));
    }

    for (int c = minCharge; c <= maxCharge; c++) {
      f_seq.push_back(charge == c ? 1.0 : 0.0); // Charge
    }
    if (enzyme_->getEnzymeType() != Enzyme::NO_ENZYME) {
      f_seq.push_back(enzyme_->isEnzymatic(peptideSeqWithFlanks.at(0), peptideSeqWithFlanks.at(2)) ? 1.0 : 0.0);
      f_seq.push_back(enzyme_->isEnzymatic(peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 3), peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 1)) ? 1.0 : 0.0);
      f_seq.push_back(static_cast<double>(enzyme_->countEnzymatic(peptideSeq)));
    }

    percolatorInNs::occurence::flankN_type flankN = peptideSeqWithFlanks.substr(0, 1);
    percolatorInNs::occurence::flankC_type flankC = peptideSeqWithFlanks.substr(peptideSeqWithFlanks.size() - 1, 1);

    // Strip peptide from termini and modifications
    std::string peptideS = peptideSeq;
    for (unsigned int ix = 0; ix < peptideSeq.size(); ++ix) {
      if (freqAA.find(peptideSeq[ix]) == string::npos) {
         peptideSeq.erase(ix--, 1);
      }
    }

    std::unique_ptr<percolatorInNs::peptideType> peptide_p(new percolatorInNs::peptideType(peptideSeq));
    // Register the ptms
    unsigned int numPTMs = 0;
    BOOST_FOREACH (const ::mzIdentML_ns::ModificationType &mod_ref, peptideMap[item.peptide_ref().get()]->Modification()) {
      BOOST_FOREACH (const ::mzIdentML_ns::CVParamType &cv_ref, mod_ref.cvParam()) {
        if (!(std::string(cv_ref.cvRef()) == "UNIMOD")) {
          ostringstream errs;
          errs << "Error: current implementation can only handle UNIMOD accessions "
               << boost::lexical_cast<string>(cv_ref.accession()) << std::endl;
          throw MyException(errs.str());
        }
        int mod_loc = boost::lexical_cast<int>(mod_ref.location());
        std::unique_ptr<percolatorInNs::modificationType> mod_p(new percolatorInNs::modificationType(mod_loc));
        if (cv_ref.accession() == "MS:1001460") {
          std::string mod_acc = "unknown";
          std::unique_ptr<percolatorInNs::freeMod> fm_p(new percolatorInNs::freeMod(mod_acc));
          mod_p->freeMod(std::move(fm_p)); // Use std::move to transfer ownership
        } else {
          int mod_acc = boost::lexical_cast<int>(cv_ref.accession().substr(7));  // Only convert text after "UNIMOD:"
          std::unique_ptr<percolatorInNs::uniMod> um_p(new percolatorInNs::uniMod(mod_acc));
          mod_p->uniMod(std::move(um_p)); // Use std::move to transfer ownership
        }
        ++numPTMs;
        peptide_p->modification().push_back(std::move(mod_p)); // Use std::move to transfer ownership
      }
    }
    
    if (po.calcPTMs) {
      f_seq.push_back(numPTMs);
    }
    if (po.pngasef) {
      f_seq.push_back(isPngasef(peptideSeqWithFlanks, isDecoy));
    }
    if (po.calcAAFrequencies) {
      computeAAFrequencies(peptideSeqWithFlanks, f_seq);
    }

    std::unique_ptr<percolatorInNs::peptideSpectrumMatch> psm_p(new percolatorInNs::peptideSpectrumMatch(
      std::move(features_p), std::move(peptide_p), psmId, isDecoy, observed_mass, theoretic_mass, charge));

    for (const std::string &proteinId : proteinIds) {
      std::unique_ptr<percolatorInNs::occurence> oc_p(new percolatorInNs::occurence(proteinId, flankN, flankC));
      psm_p->occurence().push_back(std::move(oc_p)); // Use std::move to transfer ownership
    }
    
    database->savePsm(useScanNumber, std::move(psm_p)); // Use std::move to transfer ownership
  }
  catch (std::exception const& e) {
    ostringstream temp;
    temp << "Error: parsing PSM: " << boost::lexical_cast<string>(item.id())
         << "\nThe error was: " << e.what() << std::endl;
    throw MyException(temp.str());
  }

  return;
}
