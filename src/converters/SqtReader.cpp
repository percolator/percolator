/*
    Copyright 2012 <copyright holder> <email>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#include "SqtReader.h"

SqtReader::SqtReader()
{

}

SqtReader::SqtReader(const SqtReader& other)
{

}

SqtReader::~SqtReader()
{

}

SqtReader& SqtReader::operator=(const SqtReader& other)
{
return *this;
}

bool SqtReader::operator==(const SqtReader& other) const
{
///TODO: return ...;
}

void SqtReader::readPSM(bool isDecoy, const std::string &in,  int match, const ParseOptions & po,  ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,  
			int minCharge, int maxCharge , std::string psmId , FragSpectrumScanDatabase* database ) {

  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  
  unsigned int scan;
  int charge;
  double observedMassCharge;
  double calculatedMassToCharge;

  std::istringstream instr(in), linestr;
  std::ostringstream idbuild;
  int ourPos;
  std::string line, tmp;

  double mass, deltCn, tmpdbl, otherXcorr = 0.0, xcorr = 0.0, lastXcorr =
      0.0, nSM = 0.0, tstSM = 0.0;
  bool gotL = true;
  int ms = 0;
  std::string peptide;
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::string protein;
  std::vector< std::string > proteinIds;
  //this vector is only used when we read a combined file
  std::vector< std::string > proteinIdsDecoys;
  std::map<char,int> ptmMap = po.ptmScheme; // This map should not be const declared as we use operator[], this is a FIXME for frther releases 

  while (getline(instr, line)) {
    if (line[0] == 'S') {
      linestr.clear();
      linestr.str(line);
      if (!(linestr >> tmp >> tmpdbl >> scan >> charge >> tmpdbl)) {
        cerr << "Could not parse the S line:" << endl;
        cerr << line << endl;
        exit(-1);
      }
      // Computer name might not be set, just skip this part of the line
      linestr.ignore(256, '\t');
      linestr.ignore(256, '\t');
      // First assume a MacDonald et al definition of S (9 fields)
      if (!(linestr >> observedMassCharge >> tmpdbl >> tmpdbl >> nSM)) {
        cerr << "Could not parse the S line:" << endl;
        cerr << line << endl;
        exit(-1);
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
        if (!(linestr >> tmp >> tmp >> rSp >> calculatedMassToCharge >> tmp >> xcorr >> sp
            >> matched >> expected >> peptide)) {
          cerr << "Could not parse the M line:" << endl;
          cerr << line << endl;
          exit(-1);
        }
        // replacing "*" with "-" at the terminal ends of a peptide
        if(peptide.compare(peptide.size()-1, 1, "*") == 0){
          peptide.replace(peptide.size()-1, 1, "-");
        }

        // difference between observed and calculated mass
        double dM =
            MassHandler::massDiff(observedMassCharge, calculatedMassToCharge,
                charge, peptide.substr(2, peptide.size()- 4));

        f_seq.push_back( log(max(1.0, rSp))); // rank by Sp
        f_seq.push_back( 0.0 ); // delt5Cn (leave until last M line)
        f_seq.push_back( 0.0 ); // deltCn (leave until next M line)
        f_seq.push_back( xcorr ); // Xcorr
        f_seq.push_back( sp ); // Sp
        f_seq.push_back( matched / expected ); // Fraction matched/expected ions
        f_seq.push_back( observedMassCharge ); // Observed mass
        f_seq.push_back( DataSet::peptideLength(peptide)); // Peptide length
        int nxtFeat = 8;
        for (int c = minCharge; c
        <= maxCharge; c++)
          f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge

        if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) {
          f_seq.push_back( Enzyme::isEnzymatic(peptide.at(0),peptide.at(2)) ? 1.0
              : 0.0);
          f_seq.push_back(
              Enzyme::isEnzymatic(peptide.at(peptide.size() - 3),
                  peptide.at(peptide.size() - 1))
          ? 1.0
              : 0.0);
          std::string peptid2 = peptide.substr(2, peptide.length() - 4);
          f_seq.push_back( (double)Enzyme::countEnzymatic(peptid2) );
        }
        f_seq.push_back( log(max(1.0, nSM)));
        f_seq.push_back( dM ); // obs - calc mass
        f_seq.push_back( (dM < 0 ? -dM : dM)); // abs only defined for integers on some systems
        if (po.calcPTMs) f_seq.push_back(  DataSet::cntPTMs(peptide));
        if (po.pngasef) f_seq.push_back( DataSet::isPngasef(peptide, isDecoy));
        if (po.calcAAFrequencies) {
          computeAAFrequencies(peptide, f_seq);
        }
        gotL = false;
      }
      ms++;
    }
    assert(line.size() != 0 );

    if (line.at(0) == 'L' && !gotL) {
      if (instr.peek() != 'L') gotL = true;
      
      line = line.substr(0, 2) + getRidOfUnprintables(line.substr(2));
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> protein;
      
      if(po.iscombined)
      {
	if( (line.find(po.reversedFeaturePattern) != std::string::npos) )
	  proteinIdsDecoys.push_back(protein);
	else
	  proteinIds.push_back(protein);
      }
      else 
      {
	proteinIds.push_back(protein);
      }
    }
  }

  if (xcorr > 0) {
    f_seq[1] = (xcorr - lastXcorr) / xcorr;
    f_seq[2] = (xcorr - otherXcorr) / xcorr;
  }

  if (!isfinite(f_seq[2])) std::cerr << in;

  assert(peptide.size() >= 5 );
  percolatorInNs::occurence::flankN_type flankN = peptide.substr(0,1);
  percolatorInNs::occurence::flankC_type flankC = peptide.substr(peptide.size() - 1,1);
  
  // Strip peptide from termini and modifications 
  std::string peptideSequence = peptide.substr(2, peptide.size()- 4);
  std::string peptideS = peptideSequence;
  for(unsigned int ix=0;ix<peptideSequence.size();++ix) {
    if (aaAlphabet.find(peptideSequence[ix])==string::npos && ambiguousAA.find(peptideSequence[ix])==string::npos
	&& modifiedAA.find(peptideSequence[ix])==string::npos){
      if (ptmMap.count(peptideSequence[ix])==0) {
	cerr << "Peptide sequence " << peptide << " contains modification " << peptideSequence[ix] << " that is not specified by a \"-p\" argument" << endl;
        exit(-1);
      }
      peptideSequence.erase(ix,1);
    }  
  }
  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSequence   ) );
  
  //NOTE needed for the case when a combined target/decoy file has a psm with target and decoys (copy constructor does not work)
  //std::auto_ptr< percolatorInNs::peptideType >  peptide_p2( new percolatorInNs::peptideType( peptideSequence   ) );
  
  // Register the ptms
  for(unsigned int ix=0;ix<peptideS.size();++ix) {
    if (aaAlphabet.find(peptideS[ix])==string::npos) {
      int accession = ptmMap[peptideS[ix]];
      std::auto_ptr< percolatorInNs::uniMod > um_p (new percolatorInNs::uniMod(accession));
      std::auto_ptr< percolatorInNs::modificationType >  mod_p( new percolatorInNs::modificationType(um_p,ix));
      peptide_p->modification().push_back(mod_p);      
      //peptide_p2->modification().push_back(mod_p);
      peptideS.erase(ix,1);      
    }  
  }
  
  if(!po.iscombined)
  {
    
    percolatorInNs::peptideSpectrumMatch* tmp_psm = new percolatorInNs::peptideSpectrumMatch (
	features_p,  peptide_p,psmId, isDecoy, observedMassCharge, calculatedMassToCharge, charge);
    std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(tmp_psm);

    for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) {
      std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
      psm_p->occurence().push_back(oc_p);
    }
  
    database->savePsm(scan, psm_p);
  }
  else
  { 
    std::auto_ptr< percolatorInNs::features >  features_p2( new percolatorInNs::features( *features_p->_clone() ) );
    std::auto_ptr< percolatorInNs::peptideType >  peptide_p2( new percolatorInNs::peptideType( *peptide_p->_clone() ) );
    
    if(proteinIdsDecoys.size() > 0)
    {
      percolatorInNs::peptideSpectrumMatch* tmp_psm = new percolatorInNs::peptideSpectrumMatch (
	features_p,  peptide_p,psmId, true /*is decoy*/, observedMassCharge, calculatedMassToCharge, charge);
      
      std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(tmp_psm);

      for ( std::vector< std::string >::const_iterator i = proteinIdsDecoys.begin(); i != proteinIdsDecoys.end(); ++i ) {
	std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
	psm_p->occurence().push_back(oc_p);
      }
    
      database->savePsm(scan, psm_p);
    }
    
    if(proteinIds.size() > 0)
    {      
      percolatorInNs::peptideSpectrumMatch* tmp_psm = new percolatorInNs::peptideSpectrumMatch (
	features_p2,  peptide_p2,psmId, false /*is decoy*/, observedMassCharge, calculatedMassToCharge, charge);
      
      std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(tmp_psm);

      for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) {
	std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
	psm_p->occurence().push_back(oc_p);
      }
    
      database->savePsm(scan, psm_p);
    }
    
  }
}

void SqtReader::read(const std::string fn,
    ::percolatorInNs::featureDescriptions & fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,
      bool isDecoy, const ParseOptions & po, int* maxCharge,  int* minCharge,
      parseType pType, FragSpectrumScanDatabase* database) {

  std::cerr << "Reading sqt " << fn << std::endl;
  
  std::string ptmAlphabet;
  // Normal + Amino acid + PTM + hitsPerSpectrum + doc
  const int maxNumRealFeatures = 16 + 3 + 20 * 3 + 1 + 1 + 3;

  int label;
  double *feature, *regressionFeature;
  int numSpectra;
  std::string fileId;

  int n = 0, charge = 0, ms = 0;

  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  bool look = false;
  unsigned int scanExtra;
  while (getline(sqtIn, line)) {
    if (line[0] == 'S' && sqtIn.peek() != 'S') {
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scanExtra >> charge;
      look = true;

      if ( pType == justSearchMaxMinCharge ) {
        if (*minCharge > charge) *minCharge = charge;
        if (*maxCharge < charge) *maxCharge = charge;
      }

      ms = 0;
    }
     
    if (look && line[0] == 'L' && ms < po.hitsPerSpectrum) {

	lineParse.clear();
	lineParse.str(line);
	lineParse >> tmp >> prot;
	
	if(po.iscombined)
	{
	  ++ms;
	  ++n;
	}
	else if( !isDecoy || (po.reversedFeaturePattern == "" ||
          ((line.find(po.reversedFeaturePattern, 0) != std::string::npos))))
	{
	  ++ms;
	  ++n;
	}
    }
  }
  if ( pType == justSearchMaxMinCharge ) {
    return;
  }
  if (n <= 0) {
    std::cerr << "The file " << fn << " does not contain any records"
        << std::endl;
    sqtIn.close();
    exit(-1);
  }
  sqtIn.clear();
  sqtIn.seekg(0, std::ios::beg);

  if ( fds.featureDescription().size() == 0 ) {
    addFeatureDescriptions(fds,*minCharge,
        *maxCharge,
        Enzyme::getEnzymeType() != Enzyme::NO_ENZYME 			    ,
        po.calcPTMs,
        po.pngasef,
        (po.calcAAFrequencies ? aaAlphabet : ""),
        po.calcQuadraticFeatures);
  }

  std::string seq;
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) fileId.erase(0, spos + 1);
  spos = fileId.find('.');
  if (spos != std::string::npos) fileId.erase(spos);

  std::ostringstream buff, id;
  int ix = 0, lines = 0;
  std::string scan;
  std::set<int> theMs;
  unsigned int scanNr2;

  while (getline(sqtIn, line)) {
    if (line[0] == 'S') {
      if (lines > 1) {
        readSectionS( buff.str(), fsss , theMs, isDecoy, po, *minCharge,
            *maxCharge, id.str(), database);
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
    if (line[0] == 'M') {
      ++ms;
      ++lines;
      buff << line << std::endl;
    }
    if (line[0] == 'L') {
      ++lines;
      buff << line << std::endl;
      if(po.iscombined && (int)theMs.size() < po.hitsPerSpectrum)
	theMs.insert(ms - 1);
      else if ((int)theMs.size() < po.hitsPerSpectrum &&
          ( !isDecoy || ( po.reversedFeaturePattern == "" ||
              ((line.find(po.reversedFeaturePattern, 0) != std::string::npos)))))
      {
	  theMs.insert(ms - 1);
      }
    }
  }
  if (lines > 1) {
    std::string idstr = id.str();
    readSectionS(  buff.str(), fsss, theMs, isDecoy, po,  *minCharge, *maxCharge, id.str(), database);
  }
  sqtIn.close();
}

void  SqtReader::readSectionS( std::string record , ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss, std::set<int> & theMs,  
			       bool isDecoy, const ParseOptions & po,  int minCharge, int maxCharge, std::string psmId, FragSpectrumScanDatabase* database   ) {
  std::set<int>::const_iterator it;
  for (it = theMs.begin(); it != theMs.end(); it++) {
    std::ostringstream stream;
    stream << psmId << "_" << (*it + 1);
    readPSM(isDecoy,record,*it, po, fsss, minCharge, maxCharge, stream.str(), database);
  }
  return;
}

