#include "MsgfplusReader.h"


const std::map<string, int> MsfgplusReader::msgfplusFeatures =
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

const double MsfgplusReader::neutron = 1.0033548378;  //The difference between C12 and C13

MsfgplusReader::MsfgplusReader(ParseOptions *po) : MzidentmlReader(po) {

}

MsfgplusReader::~MsfgplusReader() {



}

bool MsfgplusReader::checkValidity(const std::string &file) {
  
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



void MsfgplusReader::addFeatureDescriptions(bool doEnzyme) 
{

  push_backFeatureDescription("RawScore");
  push_backFeatureDescription("DeNovoScore");
  push_backFeatureDescription("ScoreDiff");
  push_backFeatureDescription("lnSpecEValue");
  push_backFeatureDescription("lnEValue");
  //The below are from element userParam
  push_backFeatureDescription("IsotopeError");
  push_backFeatureDescription("lnExplainedIonCurrentRatio");
  push_backFeatureDescription("lnNTermIonCurrentRatio");
  push_backFeatureDescription("lnCTermIonCurrentRatio");
  push_backFeatureDescription("lnMS2IonCurrent");
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


void MsfgplusReader::createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
        ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge,
        bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database) {

  std::auto_ptr< percolatorInNs::features > features_p(new percolatorInNs::features());
  percolatorInNs::features::feature_sequence & f_seq = features_p->feature();

  if (!item.calculatedMassToCharge().present()) {
    ostringstream temp;
    temp << "Error: calculatedMassToCharge attribute not found in PSM " 
    << boost::lexical_cast<string > (item.id())  << std::endl;
    throw MyException(temp.str());
  }

  std::string peptideSeq = peptideMap[item.peptide_ref().get()]->PeptideSequence();
  std::string peptideId = item.peptide_ref().get();
  std::vector< std::string > proteinIds;
  std::string __flankN = "";
  std::string __flankC = "";
  std::string psmid = "";

  try
  {
  
    BOOST_FOREACH(const ::mzIdentML_ns::PeptideEvidenceRefType &pepEv_ref, item.PeptideEvidenceRef()) 
    {
      std::string ref_id = pepEv_ref.peptideEvidence_ref().c_str();
      ::mzIdentML_ns::PeptideEvidenceType *pepEv = peptideEvidenceMap[ref_id];
      //NOTE check that there are not quimera peptides
      if( peptideId != std::string(pepEv->peptide_ref()))
      {
	std::cerr << "Warning : The PSM " << boost::lexical_cast<string > (item.id()) 
		  << " contains different quimera peptide sequences. "
		  << peptideMap[pepEv->peptide_ref()]->PeptideSequence() << " and " << peptideSeq 
		  << " only the proteins that contain the first peptide will be included in the PSM..\n" << std::endl;
      }
      //else
      //{
      __flankN = boost::lexical_cast<string > (pepEv->pre());
      __flankC = boost::lexical_cast<string > (pepEv->post());
      if (__flankN == "?") {__flankN = "-";} //MSGF+ sometimes outputs questionmarks here
      if (__flankC == "?") {__flankC = "-";}
      std::string proteinid = boost::lexical_cast<string > (pepEv->dBSequence_ref());
      mzIdentML_ns::SequenceCollectionType::DBSequence_type *proteinObj = proteinMap[proteinid];
      std::string proteinName = boost::lexical_cast<string > (proteinObj->accession());
      proteinIds.push_back(proteinName);
      //}
    }
    
    if(__flankC.empty() || __flankN.empty())
    {
      ostringstream temp;
      temp << "Error : The PSM " << boost::lexical_cast<string > (item.id()) << " is bad-formed." << std::endl;
      throw MyException(temp.str());
    }

    if (po->iscombined && !po->reversedFeaturePattern.empty()) {
      //NOTE taking the highest ranked PSM protein for combined search
      isDecoy = proteinIds.front().find(po->reversedFeaturePattern, 0) != std::string::npos;
    }
  
    double rank = item.rank();
    //double PI = boost::lexical_cast<double>(item.calculatedPI().get());
    int charge = item.chargeState();
    double theoretic_mass = boost::lexical_cast<double>(item.calculatedMassToCharge());
    double observed_mass = boost::lexical_cast<double>(item.experimentalMassToCharge());
    std::string peptideSeqWithFlanks = __flankN + std::string(".") + peptideSeq + std::string(".") + __flankC;
    unsigned peptide_length = peptideLength(peptideSeqWithFlanks);
    std::map<char, int> ptmMap = po->ptmScheme;
    psmid = boost::lexical_cast<string > (item.id()) + "_" + boost::lexical_cast<string > (useScanNumber) + "_" +
            boost::lexical_cast<string > (charge) + "_" + boost::lexical_cast<string > (rank);

    double RawScore = 0.0;
    double DeNovoScore = 0.0;
    double SpecEValue = 0.0;
    double EValue = 0.0;
    int IsotopeError = 0;
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
	    std::cerr << "Error  : an unmapped MS-GF+ parameter " << param_name << " was not found." << std::endl;
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
	}
    }
    
    
    //The raw theoretical mass from MSGF+ is often of the wrong isotope
    double dM = (observed_mass - (IsotopeError * neutron / charge) - theoretic_mass) / observed_mass;
    //double dM = massDiff(observed_mass, theoretic_mass, charge);  // Gives trouble because of isotopes
         //Add a small number to some logged features to avoid log(0)
    f_seq.push_back(RawScore);
    f_seq.push_back(DeNovoScore);
    f_seq.push_back(DeNovoScore - RawScore);  // Score difference (score ratio could become -inf)
    f_seq.push_back(-log(SpecEValue));
    f_seq.push_back(-log(EValue));
    f_seq.push_back(IsotopeError);
    f_seq.push_back(log(ExplainedIonCurrentRatio+0.0001));
    f_seq.push_back(log(NTermIonCurrentRatio+0.0001));
    f_seq.push_back(log(CTermIonCurrentRatio+0.0001));
    f_seq.push_back(log(MS2IonCurrent));
    f_seq.push_back(observed_mass);
    f_seq.push_back(peptideLength(peptideSeqWithFlanks));
    f_seq.push_back(dM);
    f_seq.push_back(abs(dM));
  
    for (int c = minCharge; c <= maxCharge; c++) {
      f_seq.push_back(charge == c ? 1.0 : 0.0); // Charge
    }
    if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) {
      f_seq.push_back(Enzyme::isEnzymatic(peptideSeqWithFlanks.at(0), peptideSeqWithFlanks.at(2)) ? 1.0 : 0.0);
      f_seq.push_back(Enzyme::isEnzymatic(peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 3), peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 1)) ? 1.0 : 0.0);
      f_seq.push_back((double) Enzyme::countEnzymatic(peptideSeq));
    }

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
	   ostringstream temp;
          temp << "Error : Peptide sequence " << peptideSeqWithFlanks
                  << " contains modification " << peptideSeq[ix] << " that is not specified by a \"-p\" argument" << std::endl;
          throw MyException(temp.str());
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
  }
  // Try-Catch statement to find potential errors among the features.
  catch(std::exception const& e)
  {
    ostringstream temp;
    temp << "Error : parsing PSM: " << boost::lexical_cast<string > (item.id()) 
    << "\nThe error was: " << e.what() << std::endl;
    throw MyException(temp.str());
  }
 
  return;
}
