#include "SqtReader.h"

SqtReader::SqtReader(ParseOptions *po):Reader(po)
{
}

SqtReader::~SqtReader()
{

}

void SqtReader::readPSM(bool isDecoy, const std::string &in,int match,  
			std::string psmId,boost::shared_ptr<FragSpectrumScanDatabase> database) 
{

  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  unsigned int scan;
  int charge;
  double observedMassCharge;
  double calculatedMassToCharge;

  std::istringstream instr(in), linestr;
  std::ostringstream idbuild;
  int ourPos;
  std::string line, tmp;

  double mass, deltCn, tmpdbl, otherXcorr = 0.0, xcorr = 0.0, lastXcorr = 0.0, nSM = 0.0, tstSM = 0.0;
  bool gotL = true;
  int ms = 0;
  std::string peptide;
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::string protein;
  std::vector< std::string > proteinIds;
  std::map<char,int> ptmMap = po->ptmScheme; 

  while (getline(instr, line)) 
  {
    if (line[0] == 'S') 
    {
      linestr.clear();
      linestr.str(line);
      if (!(linestr >> tmp >> tmpdbl >> scan >> charge >> tmpdbl)) 
      {
	ostringstream temp;
        temp << "Error : can not parse the S line: " << line << endl;
	throw MyException(temp.str());
      }
      // Computer name might not be set, just skip this part of the line
      linestr.ignore(256, '\t');
      linestr.ignore(256, '\t');
      // First assume a MacDonald et al definition of S (9 fields)
      if (!(linestr >> observedMassCharge >> tmpdbl >> tmpdbl >> nSM)) 
      {
	ostringstream temp;
        temp << "Error : can not parse the S line: " << line << endl;
	throw MyException(temp.str());
      }
      // Check if the Yate's lab definition (10 fields) is valid
      // http://fields.scripps.edu/sequest/SQTFormat.html
      //
      if (linestr >> tstSM) 
      {
        nSM = tstSM;
      }
    }
    if (line[0] == 'M') 
    {
      linestr.clear();
      linestr.str(line);

      if (ms == 1) 
      {
        linestr >> tmp >> tmp >> tmp >> tmp >> deltCn >> otherXcorr;
        lastXcorr = otherXcorr;
      } 
      else 
      {
        linestr >> tmp >> tmp >> tmp >> tmp >> tmp >> lastXcorr;
      }
      if (match == ms) 
      {
        double rSp, sp, matched, expected;
        linestr.seekg(0, ios::beg);
        if (!(linestr >> tmp >> tmp >> rSp >> calculatedMassToCharge >> tmp >> xcorr >> sp >> matched >> expected >> peptide)) 
	{
	  ostringstream temp;
	  temp << "Error : can not parse the M line: " << line << endl;
	  throw MyException(temp.str());
        }
        
        // difference between observed and calculated mass
        double dM = massDiff(observedMassCharge, calculatedMassToCharge,charge);
	
        f_seq.push_back( log(max(1.0, rSp))); // rank by Sp
        f_seq.push_back( 0.0 ); // delt5Cn (leave until last M line)
        f_seq.push_back( 0.0 ); // deltCn (leave until next M line)
        f_seq.push_back( xcorr ); // Xcorr
        f_seq.push_back( sp ); // Sp
        f_seq.push_back( matched / expected ); // Fraction matched/expected ions
        f_seq.push_back( observedMassCharge ); // Observed mass
        f_seq.push_back(peptideLength(peptide)); // Peptide length
        int nxtFeat = 8;
        for (int c = minCharge; c <= maxCharge; c++)
          f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge

        if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) 
	{
          f_seq.push_back( Enzyme::isEnzymatic(peptide.at(0),peptide.at(2)) ? 1.0 : 0.0);
          f_seq.push_back( Enzyme::isEnzymatic(peptide.at(peptide.size() - 3),peptide.at(peptide.size() - 1)) ? 1.0 : 0.0);
          std::string peptid2 = peptide.substr(2, peptide.length() - 4);
          f_seq.push_back( (double)Enzyme::countEnzymatic(peptid2) );
        }
        f_seq.push_back( log(max(1.0, nSM)));
        f_seq.push_back( dM ); // obs - calc mass
        f_seq.push_back( abs(dM) ); // abs only defined for integers on some systems
        if (po->calcPTMs) 
	{
	  f_seq.push_back(cntPTMs(peptide));
	}
	if (po->pngasef) 
	{
	  f_seq.push_back(isPngasef(peptide, isDecoy));
	}
        if (po->calcAAFrequencies) 
	{
          computeAAFrequencies(peptide, f_seq);
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

  if (xcorr > 0) {
    f_seq[1] = (xcorr - lastXcorr) / xcorr;
    f_seq[2] = (xcorr - otherXcorr) / xcorr;
  }
  
  percolatorInNs::occurence::flankN_type flankN = peptide.substr(0,1);
  percolatorInNs::occurence::flankC_type flankC = peptide.substr(peptide.size() - 1,1);
  
  // Strip peptide from termini and modifications 
  std::string peptideSequence = peptide.substr(2, peptide.size()- 4);
  std::string peptideS = peptideSequence;
  for(unsigned int ix=0;ix<peptideSequence.size();++ix) 
  {
    if (aaAlphabet.find(peptideSequence[ix])==string::npos && 
	ambiguousAA.find(peptideSequence[ix])==string::npos && 
	additionalAA.find(peptideSequence[ix])==string::npos)
    {
      if (ptmMap.count(peptideSequence[ix])==0) 
      {
	ostringstream temp;
	temp << "Error : Peptide sequence " << peptide << " contains modification " 
	<< peptideSequence[ix] << " that is not specified by a \"-p\" argument" << endl;
        throw MyException(temp.str());
      }
      peptideSequence.erase(ix,1);
    }  
  }
  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSequence   ) );
  // Register the ptms
  for(unsigned int ix=0;ix<peptideS.size();++ix) 
  {
    if (aaAlphabet.find(peptideS[ix])==string::npos && 
	ambiguousAA.find(peptideS[ix])==string::npos && 
	additionalAA.find(peptideS[ix])==string::npos)
    {
      int accession = ptmMap[peptideS[ix]];
      std::auto_ptr< percolatorInNs::uniMod > um_p (new percolatorInNs::uniMod(accession));
      std::auto_ptr< percolatorInNs::modificationType >  mod_p( new percolatorInNs::modificationType(um_p,ix));
      peptide_p->modification().push_back(mod_p);      
      peptideS.erase(ix,1);      
    }  
  }
  
  if(po->iscombined)
  {
    //NOTE when combine search the PSM will take the identity of its first ranked protein
    isDecoy = proteinIds.front().find(po->reversedFeaturePattern, 0) != std::string::npos;
  }
  
  std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  
  psm_p(new percolatorInNs::peptideSpectrumMatch
  (features_p,  peptide_p,psmId, isDecoy, observedMassCharge, calculatedMassToCharge, charge));

  for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) 
  {
    std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
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
      minCharge = std::min(minCharge,charge);
      maxCharge = std::max(maxCharge,charge);
    }
     
  }
  sqtIn.close();
}


void SqtReader::read(const std::string &fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database) 
{

  std::string ptmAlphabet;
  int label;
  double *feature, *regressionFeature;
  int numSpectra;
  std::string fileId;
  int charge;
  int n = 0, ms = 0;
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) 
  {
    ostringstream temp;
    temp << "Error : can not open file " << fn << std::endl;
    throw MyException(temp.str());
  }

  std::string seq;
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) 
  {  
    fileId.erase(0, spos + 1);
  }
  spos = fileId.find('.');
  if (spos != std::string::npos)
  {
    fileId.erase(spos);
  }
  std::ostringstream buff, id;
  int ix = 0, lines = 0;
  std::string scan;
  std::set<int> theMs;
  unsigned int scanNr2;

  while (getline(sqtIn, line)) 
  {
    if (line[0] == 'S') 
    {
      if (lines > 1) 
      {
        readSectionS( buff.str(), theMs, isDecoy,id.str(), database);
      }
      buff.str("");
      buff.clear();
      id.str("");
      lines = 1;
      buff << line << std::endl;
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scan >> charge;
      id << fileId << '_' << scan << '_' << charge;
      ms = 0;
      theMs.clear();
    }
    if (line[0] == 'M') 
    {
      ++ms;
      ++lines;
      buff << line << std::endl;
    }
    if (line[0] == 'L') 
    {
      ++lines;
      buff << line << std::endl;
     if ((int)theMs.size() < po->hitsPerSpectrum && 
       ( !isDecoy || ( po->reversedFeaturePattern == "" || 
       ((line.find(po->reversedFeaturePattern, 0) != std::string::npos)))))
      {
	  theMs.insert(ms - 1);
      }
    }
  }
  if (lines > 1) 
  {
    readSectionS(buff.str(), theMs, isDecoy, id.str(), database);
  }
  sqtIn.close();
}

void  SqtReader::readSectionS(const std::string &record,std::set<int> & theMs, bool isDecoy,
			       std::string psmId,boost::shared_ptr<FragSpectrumScanDatabase> database) 
{
  std::set<int>::const_iterator it;
  for (it = theMs.begin(); it != theMs.end(); it++) 
  {
    std::ostringstream stream;
    stream << psmId << "_" << (*it + 1);
    readPSM(isDecoy,record,*it, stream.str(), database);
  }
  return;
}

bool SqtReader::checkValidity(const std::string &file)
{
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
  push_backFeatureDescription("lnrSp");
  push_backFeatureDescription("deltLCn");
  push_backFeatureDescription("deltCn");
  push_backFeatureDescription("Xcorr");
  push_backFeatureDescription("Sp");
  push_backFeatureDescription("IonFrac");
  push_backFeatureDescription("Mass");
  push_backFeatureDescription("PepLen");

  for (int charge = minCharge; charge <= maxCharge; charge++) 
  {
    std::ostringstream cname;
    cname << "Charge" << charge;
    push_backFeatureDescription(cname.str().c_str());

  }
  if (doEnzyme) 
  {
    push_backFeatureDescription("enzN");
    push_backFeatureDescription("enzC");
    push_backFeatureDescription("enzInt");
  }
  push_backFeatureDescription("lnNumSP");
  push_backFeatureDescription("dM");
  push_backFeatureDescription("absdM");
  if (po->calcPTMs) 
  {
    push_backFeatureDescription("ptm");
  }
  if (po->pngasef) 
  {
    push_backFeatureDescription("PNGaseF");
  }
  if (po->calcAAFrequencies)
  {
    for (std::string::const_iterator it = aaAlphabet.begin(); it != aaAlphabet.end(); it++)
    {
      std::string temp = boost::lexical_cast<std::string>(*it)+"-Freq";
      push_backFeatureDescription(temp.c_str());
    }
  }
}