#include "SqtReader.h"

//default score vector //TODO move this to a file or input parameter                          
const std::map<string,double> SqtReader::sqtFeaturesDefaultValue =
boost::assign::map_list_of("lnrSp", 0.0)
("deltLCn",0.0)
("deltCn", 1.61)
("Xcorr", 1.1)
("Sp", 0.0)
("IonFrac", 0.0)
("Mass", 0.0)
("PepLen", -0.573)
("Charge1", 0.0335)
("Charge2", 0.149)
("Charge3", -0.156);

SqtReader::SqtReader(ParseOptions po):Reader(po)
{
}

SqtReader::~SqtReader()
{

}

void SqtReader::readPSM(bool isDecoy, const std::string &in, int match,  
			std::string& fileId, boost::shared_ptr<FragSpectrumScanDatabase> database) {
  std::unique_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  unsigned int scan;
  int charge;
  double observedMassCharge;
  double calculatedMassToCharge;

  std::istringstream instr(in), linestr;
  std::ostringstream idbuild;
  std::string line, tmp;

  double deltCn, tmpdbl, otherXcorr = 0.0, xcorr = 0.0, lastXcorr = 0.0, nSM = 0.0, tstSM = 0.0;
  bool gotL = true;
  int ms = 0;
  std::string peptide, peptideNoMods;
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::string protein;
  std::vector< std::string > proteinIds;
  std::map<char,int> ptmMap = po.ptmScheme; 

  while (getline(instr, line)) {
    if (line[0] == 'S') {
      linestr.clear();
      linestr.str(line);
      if (!(linestr >> tmp >> tmpdbl >> scan >> charge >> tmpdbl)) {
	      ostringstream temp;
        temp << "Error : can not parse the S line: " << line << endl;
	      throw MyException(temp.str());
      }
      // Computer name might not be set, just skip this part of the line
      linestr.ignore(256, '\t');
      linestr.ignore(256, '\t');
      // First assume a MacDonald et al definition of S (9 fields)
      if (!(linestr >> observedMassCharge >> tmpdbl >> tmpdbl >> nSM)) {
      	ostringstream temp;
        temp << "Error : can not parse the S line: " << line << endl;
	      throw MyException(temp.str());
      }
      // Check if the Yate's lab definition (10 fields) is valid
      // http://fields.scripps.edu/sequest/SQTFormat.html
      //
      if (linestr >> tstSM) {
        nSM = tstSM;
      }
    }
    if (line[0] == 'M') {
      linestr.clear();
      linestr.str(line);

      if (ms == 1) {
        linestr >> tmp >> tmp >> tmp >> tmp >> deltCn >> otherXcorr;
        lastXcorr = otherXcorr;
      } else {
        linestr >> tmp >> tmp >> tmp >> tmp >> tmp >> lastXcorr;
      }
      if (match == ms) {
        double rSp, sp, matched, expected;
        linestr.seekg(0, ios::beg);
        if (!(linestr >> tmp >> tmp >> rSp >> calculatedMassToCharge >> tmp >> 
                         xcorr >> sp >> matched >> expected >> peptide)) {
	        ostringstream temp;
	        temp << "Error : can not parse the M line: " << line << endl;
	        throw MyException(temp.str());
        }
        
        peptideNoMods = removePTMs(peptide, ptmMap);
        // difference between observed and calculated mass
        double dM = massDiff(observedMassCharge, 
          calculatedMassToCharge,static_cast<unsigned int>(charge));
	
        f_seq.push_back( log(max(1.0, rSp))); // rank by Sp
        f_seq.push_back( 0.0 ); // delt5Cn (leave until last M line)
        f_seq.push_back( 0.0 ); // deltCn (leave until next M line)
        f_seq.push_back( xcorr ); // Xcorr
        f_seq.push_back( sp ); // Sp
        f_seq.push_back( matched / expected ); // Fraction matched/expected ions
        f_seq.push_back( observedMassCharge ); // Observed mass
        f_seq.push_back(peptideLength(peptide)); // Peptide length
        for (int c = minCharge; c <= maxCharge; c++)
          f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge

        if (enzyme_->getEnzymeType() != Enzyme::NO_ENZYME) {
          f_seq.push_back( enzyme_->isEnzymatic(peptideNoMods.at(0),peptideNoMods.at(2)) ? 1.0 : 0.0);
          f_seq.push_back( enzyme_->isEnzymatic(peptideNoMods.at(peptideNoMods.size() - 3),peptideNoMods.at(peptideNoMods.size() - 1)) ? 1.0 : 0.0);
          std::string peptide2 = peptideNoMods.substr(2, peptideNoMods.length() - 4);
          f_seq.push_back( (double)enzyme_->countEnzymatic(peptide2) );
        }
        
        f_seq.push_back( log(max(1.0, nSM)));
        f_seq.push_back( dM ); // obs - calc mass
        f_seq.push_back( abs(dM) ); // abs only defined for integers on some systems
        
        if (po.calcPTMs) {
	        f_seq.push_back(cntPTMs(peptide, ptmMap));
	      }
	      if (po.pngasef) {
	        f_seq.push_back(isPngasef(peptide, isDecoy));
	      }
        if (po.calcAAFrequencies) {
          computeAAFrequencies(peptideNoMods, f_seq);
        }
        gotL = false;
      }
      ms++;
    }
    assert(line.size() != 0 );

    if (line.at(0) == 'L' && !gotL) 
    {
      if (instr.peek() != 'L') gotL = true;
      
      line = line.substr(0, 2) + getRidOfUnprintables(line.substr(2));
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> protein;
      
      proteinIds.push_back(protein);
      
    }
  }
  f_seq[1] = (xcorr - lastXcorr) / (std::max)(1.0,xcorr); // delt5Cn
  f_seq[2] = (xcorr - otherXcorr) / (std::max)(1.0,xcorr); // deltCn
  
  percolatorInNs::occurence::flankN_type flankN = peptide.substr(0,1);
  percolatorInNs::occurence::flankC_type flankC = peptide.substr(peptide.size() - 1,1);
  
  // Strip peptide from termini and modifications 
  std::string peptideSequence = peptide.substr(2, peptide.size()- 4);
  std::string peptideSeqNoMods = peptideNoMods.substr(2, peptideNoMods.size()- 4);
  std::unique_ptr< percolatorInNs::peptideType > peptide_p( new percolatorInNs::peptideType(peptideSeqNoMods) );
  // Register the ptms
  for (unsigned int ix = 0;ix < peptideSequence.size();++ix) {
    if (freqAA.find(peptideSequence[ix]) == string::npos) {
      std::unique_ptr< percolatorInNs::modificationType >  mod_p( new percolatorInNs::modificationType(static_cast<int>(ix)));
      if (peptideSequence[ix] == '[') {
        unsigned int posEnd = static_cast<unsigned int>(peptideSequence.substr(ix).find_first_of(']'));
        std::string modAcc = peptideSequence.substr(ix + 1, posEnd - 1);
        std::unique_ptr< percolatorInNs::freeMod > fm_p (new percolatorInNs::freeMod(modAcc));
        mod_p->freeMod(fm_p);
        peptideSequence.erase(ix--, posEnd + 1);
      } else {
        int accession = ptmMap[peptideSequence[ix]];
        std::unique_ptr< percolatorInNs::uniMod > um_p (new percolatorInNs::uniMod(accession));
        mod_p->uniMod(um_p);  
        peptideSequence.erase(ix--,1);
      }
      peptide_p->modification().push_back(mod_p);    
    }  
  }
  
  if (po.iscombined) {
    // if one of the proteins is a target protein, we classify the PSM as a target PSM
    isDecoy = true;
    for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) {
      if (i->find(po.reversedFeaturePattern, 0) == std::string::npos) {
        isDecoy = false;
        break;
      }
    }
  }
  
  unsigned int rank = static_cast<unsigned int>(match + 1);
  std::string psmId = createPsmId(fileId, observedMassCharge, scan, charge, rank);
  
  std::unique_ptr< percolatorInNs::peptideSpectrumMatch > psm_p(
      new percolatorInNs::peptideSpectrumMatch(features_p, peptide_p, psmId, 
          isDecoy, observedMassCharge, calculatedMassToCharge, charge));

  for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) {
    std::unique_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
    psm_p->occurence().push_back(oc_p);
  }
  
  database->savePsm(scan, psm_p);
}

void SqtReader::getMaxMinCharge(const std::string &fn, bool isDecoy)
{
  int charge = 0;
  std::string line;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  unsigned int scanExtra;
  std::string tmp;
  
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) 
  {
    ostringstream temp;
    temp << "Error : can not open file " << fn << std::endl;
    throw MyException(temp.str());
  }

  while (getline(sqtIn, line)) 
  {
    if (line[0] == 'S' && sqtIn.peek() != 'S') 
    {
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scanExtra >> charge;
      minCharge = (std::min)(minCharge,charge);
      maxCharge = (std::max)(maxCharge,charge);
    }
     
  }
  sqtIn.close();
}


void SqtReader::read(const std::string &fn, bool isDecoy, 
    boost::shared_ptr<FragSpectrumScanDatabase> database) {
  std::string ptmAlphabet;
  std::string fileId;
  int ms = 0;
  std::string line, tmp, prot;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) {
    ostringstream temp;
    temp << "Error : can not open file " << fn << std::endl;
    throw MyException(temp.str());
  }

  std::string seq;
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) {
    fileId.erase(0, spos + 1);
  }
  spos = fileId.find('.');
  if (spos != std::string::npos) {
    fileId.erase(spos);
  }
  std::ostringstream buff;
  int lines = 0;
  std::set<int> theMs;

  while (getline(sqtIn, line)) {
    if (line[0] == 'S') {
      if (lines > 1) {
        readSectionS( buff.str(), theMs, isDecoy, fileId, database);
      }
      buff.str("");
      buff.clear();
      lines = 1;
      buff << line << std::endl;
      ms = 0;
      theMs.clear();
    }
    if (line[0] == 'M') {
      ++ms;
      ++lines;
      buff << line << std::endl;
    }
    if (line[0] == 'L') {
      ++lines;
      buff << line << std::endl;
      if ((int)theMs.size() < po.hitsPerSpectrum && 
           ( !po.iscombined || !isDecoy || ( po.reversedFeaturePattern == "" || 
           ((line.find(po.reversedFeaturePattern, 0) != std::string::npos))))) {
	      theMs.insert(ms - 1);
      }
    }
  }
  if (lines > 1) {
    readSectionS(buff.str(), theMs, isDecoy, fileId, database);
  }
  sqtIn.close();
}

void SqtReader::readSectionS(const std::string &record, std::set<int>& theMs, bool isDecoy,
			       std::string& fileId, boost::shared_ptr<FragSpectrumScanDatabase> database) {
  std::set<int>::const_iterator it;
  for (it = theMs.begin(); it != theMs.end(); it++) {
    readPSM(isDecoy, record, *it, fileId, database);
  }
  return;
}

bool SqtReader::checkValidity(const std::string &file) {
  bool isvalid = true;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  if (!fileIn) 
  {
    ostringstream temp;
    temp << "Error : can not open file " << file << std::endl;
    throw MyException(temp.str());
  }
  std::string line;
  if (!getline(fileIn, line)) 
  {
    ostringstream temp;
    temp << "Error : can not read file " << file << std::endl;
    throw MyException(temp.str());
  }
  fileIn.close();
  if (line.find("SQTGenerator") == std::string::npos) 
  {
    ostringstream temp;
    temp << "Error : SQT file not generated by CRUX: " << file << std::endl;
    throw MyException(temp.str());
  }
 
  return isvalid;
}

bool SqtReader::checkIsMeta(const std::string &file)
{
  //NOTE assuming the file has been tested before
  bool isMeta;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  std::string line;
  getline(fileIn, line);
  fileIn.close();
  if (line.size() > 1 && line[0]=='H' && (line[1]=='\t' || line[1]==' '))
  {
    isMeta = false;
  }
  else
  {
    isMeta = true;
  }
  return isMeta;
}

void SqtReader::addFeatureDescriptions(bool doEnzyme) 
{
  push_backFeatureDescription("lnrSp","",sqtFeaturesDefaultValue.at("lnrSp"));
  push_backFeatureDescription("deltLCn","",sqtFeaturesDefaultValue.at("deltLCn"));
  push_backFeatureDescription("deltCn","",sqtFeaturesDefaultValue.at("deltCn"));
  push_backFeatureDescription("Xcorr","",sqtFeaturesDefaultValue.at("Xcorr"));
  push_backFeatureDescription("Sp","",sqtFeaturesDefaultValue.at("Sp"));
  push_backFeatureDescription("IonFrac","",sqtFeaturesDefaultValue.at("IonFrac"));
  push_backFeatureDescription("Mass","",sqtFeaturesDefaultValue.at("Mass"));
  push_backFeatureDescription("PepLen","",sqtFeaturesDefaultValue.at("PepLen"));

  for (int charge = minCharge; charge <= maxCharge; charge++) 
  {
    std::ostringstream cname;
    cname << "Charge" << charge;
    double value = 0.0;
    //UGLY!!
    if(charge == 1 || charge == 2 || charge == 3)  {
         value = sqtFeaturesDefaultValue.at(cname.str());
    }
    push_backFeatureDescription(cname.str().c_str(),"",value);

  }
  
  //the rest of the values will have a default of 0.0
  if (doEnzyme) 
  {
    push_backFeatureDescription("enzN");
    push_backFeatureDescription("enzC");
    push_backFeatureDescription("enzInt");
  }
  push_backFeatureDescription("lnNumSP");
  push_backFeatureDescription("dM");
  push_backFeatureDescription("absdM");
  if (po.calcPTMs) 
  {
    push_backFeatureDescription("ptm");
  }
  if (po.pngasef) 
  {
    push_backFeatureDescription("PNGaseF");
  }
  if (po.calcAAFrequencies)
  {
    BOOST_FOREACH (const char aa, freqAA) {
      std::string temp = std::string(1,aa) + "-Freq";
      push_backFeatureDescription(temp.c_str());
    }
  }
}
