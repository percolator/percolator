#include "SqtReader.h"

SqtReader::SqtReader(ParseOptions po):Reader(po)
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
  std::map<char,int> ptmMap = po.ptmScheme; 

  while (getline(instr, line)) 
  {
    if (line[0] == 'S') 
    {
      linestr.clear();
      linestr.str(line);
      if (!(linestr >> tmp >> tmpdbl >> scan >> charge >> tmpdbl)) 
      {
        cerr << "Could not parse the S line:" << endl;
        cerr << line << endl;
        exit(-1);
      }
      // Computer name might not be set, just skip this part of the line
      linestr.ignore(256, '\t');
      linestr.ignore(256, '\t');
      // First assume a MacDonald et al definition of S (9 fields)
      if (!(linestr >> observedMassCharge >> tmpdbl >> tmpdbl >> nSM)) 
      {
        cerr << "Could not parse the S line:" << endl;
        cerr << line << endl;
        exit(-1);
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
          cerr << "Could not parse the M line:" << endl;
          cerr << line << endl;
          exit(-1);
        }
        // replacing "*" with "-" at the terminal ends of a peptide
        if(peptide.compare(peptide.size()-1, 1, "*") == 0)
	{
          peptide.replace(peptide.size()-1, 1, "-");
        }
        // difference between observed and calculated mass
        double dM = MassHandler::massDiff(observedMassCharge, calculatedMassToCharge,charge, peptide.substr(2, peptide.size()- 4));

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
        if (po.calcPTMs) 
	{
	  f_seq.push_back(cntPTMs(peptide));
	}
	if (po.pngasef) 
	{
	  f_seq.push_back(isPngasef(peptide, isDecoy));
	}
        if (po.calcAAFrequencies) 
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

  if (!isfinite(f_seq[2])) std::cerr << in;

  assert(peptide.size() >= po.peptidelength );
  percolatorInNs::occurence::flankN_type flankN = peptide.substr(0,1);
  percolatorInNs::occurence::flankC_type flankC = peptide.substr(peptide.size() - 1,1);
  
  // Strip peptide from termini and modifications 
  std::string peptideSequence = peptide.substr(2, peptide.size()- 4);
  std::string peptideS = peptideSequence;
  for(unsigned int ix=0;ix<peptideSequence.size();++ix) 
  {
    if (aaAlphabet.find(peptideSequence[ix])==string::npos && ambiguousAA.find(peptideSequence[ix])==string::npos && modifiedAA.find(peptideSequence[ix])==string::npos)
    {
      if (ptmMap.count(peptideSequence[ix])==0) 
      {
	cerr << "Peptide sequence " << peptide << " contains modification " << peptideSequence[ix] << " that is not specified by a \"-p\" argument" << endl;
        exit(-1);
      }
      peptideSequence.erase(ix,1);
    }  
  }
  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSequence   ) );
  // Register the ptms
  for(unsigned int ix=0;ix<peptideS.size();++ix) 
  {
    if (aaAlphabet.find(peptideS[ix])==string::npos) 
    {
      int accession = ptmMap[peptideS[ix]];
      std::auto_ptr< percolatorInNs::uniMod > um_p (new percolatorInNs::uniMod(accession));
      std::auto_ptr< percolatorInNs::modificationType >  mod_p( new percolatorInNs::modificationType(um_p,ix));
      peptide_p->modification().push_back(mod_p);      
      peptideS.erase(ix,1);      
    }  
  }
  
  if(po.iscombined)
  {
    //NOTE when combine search the PSM will take the identity of its first ranked Peptide
    isDecoy = proteinIds.front().find(po.reversedFeaturePattern, 0) != std::string::npos;
  }
  
  percolatorInNs::peptideSpectrumMatch* tmp_psm = new percolatorInNs::peptideSpectrumMatch (
	features_p,  peptide_p,psmId, isDecoy, observedMassCharge, calculatedMassToCharge, charge);
  std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(tmp_psm);

  for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) 
  {
    std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
    psm_p->occurence().push_back(oc_p);
  }
  
  database->savePsm(scan, psm_p);
  
}

void SqtReader::getMaxMinCharge(string fn)
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
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }

  while (getline(sqtIn, line)) 
  {
    if (line[0] == 'S' && sqtIn.peek() != 'S') 
    {
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scanExtra >> charge;
      if (minCharge > charge) minCharge = charge;
      if (maxCharge < charge) maxCharge = charge; 
    }
     
  }
  sqtIn.close();
}


void SqtReader::read(const std::string fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database) 
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
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  
  bool look = false;
  unsigned int scanExtra;
  while (getline(sqtIn, line)) 
  {
    if (line[0] == 'S' && sqtIn.peek() != 'S') 
    {
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scanExtra >> charge;
      look = true;
      ms = 0;
    }
     
    if (look && line[0] == 'L' && ms < po.hitsPerSpectrum) 
    {
	lineParse.clear();
	lineParse.str(line);
	lineParse >> tmp >> prot;
	if( !isDecoy || (po.reversedFeaturePattern == "" || ( (line.find(po.reversedFeaturePattern, 0) != std::string::npos) ) ) )
	{
	  ++ms;
	  ++n;
	}
    }
  }

  if (n <= 0) 
  {
    std::cerr << "The file " << fn << " does not contain any records"
        << std::endl;
    sqtIn.close();
    exit(-1);
  }
  sqtIn.clear();
  sqtIn.seekg(0, std::ios::beg);

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
     if ((int)theMs.size() < po.hitsPerSpectrum && ( !isDecoy || ( po.reversedFeaturePattern == "" || ((line.find(po.reversedFeaturePattern, 0) != std::string::npos)))))
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

void  SqtReader::readSectionS( std::string record,std::set<int> & theMs, bool isDecoy,
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

bool SqtReader::checkValidity(const std::string file)
{
  bool ismeta = true;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  if (!fileIn) 
  {
    std::cerr << "Could not open file " << file << std::endl;
    exit(-1);
  }
  std::string line;
  if (!getline(fileIn, line)) 
  {
    std::cerr << "Could not read file " << file << std::endl;
    exit(-1);
  }
  fileIn.close();
  if (line.size() > 1 && line[0]=='H' && (line[1]=='\t' || line[1]==' '))
  {
    if (line.find("SQTGenerator") == std::string::npos) 
    {
      std::cerr << "SQT file not generated by CRUX: " << file << std::endl;
      exit(-1);
    }
  }
  else
  {
    ismeta = false;
  }
  
  return ismeta;
}

void SqtReader::addFeatureDescriptions(bool doEnzyme,const std::string& aaAlphabet,std::string fn) 
{
  push_backFeatureDescription("lnrSp");
  push_backFeatureDescription("deltLCn");
  push_backFeatureDescription("deltCn");
  push_backFeatureDescription("Xcorr");
  push_backFeatureDescription("Sp");
  push_backFeatureDescription("IonFrac");
  push_backFeatureDescription("Mass");
  push_backFeatureDescription("PepLen");

  for (int charge = minCharge; charge <= maxCharge; ++charge) 
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
    for (std::string::const_iterator it = aaAlphabet.begin(); it != aaAlphabet.end(); it++)
    {
      std::string temp = boost::lexical_cast<std::string>(*it)+"-Freq";
      push_backFeatureDescription(temp.c_str());
    }
  }
  
  //NOTE this is not being filled up later in the PSM
  /**if (po.calcQuadraticFeatures) 
  {
    for (int f1 = 1; f1 < f_seq.featureDescription().size(); ++f1) 
    {
      for (int f2 = 0; f2 < f1; ++f2) 
      {
        std::ostringstream feat;
        feat << "f" << f1 + 1 << "*" << "f" << f2 + 1;
        push_backFeatureDescription(feat.str().c_str());
      }
    }
  }**/
}

void SqtReader::readRetentionTime(string filename) 
{
  MSReader r;
  Spectrum s;
  r.setFilter(MS2);
  char* cstr = new char[filename.size() + 1];
  strcpy(cstr, filename.c_str());
  // read first spectrum
  r.readFile(cstr, s);
  while(s.getScanNumber() != 0)
  {
    // check whether an EZ lines is available
    if(s.sizeEZ() != 0)
    {
      // for each EZ line (each psm)
      for(int i = 0; i<s.sizeEZ(); i++)
      {
        // save experimental mass and retention time
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).mh);
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).pRTime);
      }
    }
    // if no EZ line is available, check for an RTime lines
    else if((double)s.getRTime() != 0)
    {
      scan2rt[s.getScanNumber()].push_back(s.getRTime());
    }
    // if neither EZ nor I lines are available
    else
    {
      cout << "The ms2 in input does not appear to contain retention time "
          << "information. Please run without -2 option.";
      exit(-1);
    }
    // read next scan
    r.readFile(NULL, s);
  }
  delete[] cstr;
}

void SqtReader::storeRetentionTime(boost::shared_ptr<FragSpectrumScanDatabase> database)
{
  // for each spectra from the ms2 file
  typedef std::map<int, vector<double> > map_t;
  BOOST_FOREACH(map_t::value_type& i, scan2rt)
  {
    // scan number
    int scanNr = i.first;
    // related retention times
    vector<double>* rTimes = &(i.second);
    if(database->getFSS(scanNr).get()!=0)
    {
      fragSpectrumScan fss = *(database->getFSS(scanNr));
      fragSpectrumScan::peptideSpectrumMatch_sequence& psmSeq = fss.peptideSpectrumMatch();
      // retention time to be stored
      double storeMe = 0;
      // if rTimes only contains one element
      if(rTimes->size()==1)
      {
        // take that as retention time
        storeMe = rTimes->at(0);
      }
      else
      {
        // else, take retention time of psm that has observed mass closest to
        // theoretical mass (smallest massDiff)
        double massDiff = std::numeric_limits<double>::max(); // + infinity
        for (fragSpectrumScan::peptideSpectrumMatch_iterator psmIter_i = psmSeq.begin(); psmIter_i != psmSeq.end(); ++psmIter_i) 
	{
          // skip decoy
          if(psmIter_i->isDecoy() != true)
	  {
            double cm = psmIter_i->calculatedMassToCharge();
            double em = psmIter_i->experimentalMassToCharge();
            // if a psm with observed mass closer to theoretical mass is found
            if(abs(cm-em) < massDiff)
	    {
              // update massDiff
              massDiff = abs(cm-em);
              // get corresponding retention time
              vector<double>::const_iterator r = rTimes->begin();
              for(; r<rTimes->end(); r=r+2)
	      {
                double rrr = *r;
                double exm = psmIter_i->experimentalMassToCharge();
                if(*r==psmIter_i->experimentalMassToCharge())
		{
                  storeMe = *(r+1);
                  r = rTimes->end();
                }
              }
            }
          }
        }
      }
      // store retention time for all psms in fss
      for (fragSpectrumScan::peptideSpectrumMatch_iterator psmIter = psmSeq.begin(); psmIter != psmSeq.end(); ++psmIter) 
      {
        psmIter->observedTime().set(storeMe);
      }
      database->putFSS(fss);
    }
  }
}
