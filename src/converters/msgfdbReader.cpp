#include "msgfdbReader.h"

msgfdbReader::msgfdbReader(ParseOptions po):Reader(po)
{
  
}

msgfdbReader::~msgfdbReader()
{

}


//Some functions to handle strings
std::vector<std::string> &msgfdbReader::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
      if(item.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")!=std::string::npos){
	elems.push_back(item); 
      }
    }
    return elems;
}

std::vector<std::string> msgfdbReader::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

/**
void msgfdbReader::remove_endl(std::string s)
{
  int pos=s.find('\n');
  if(pos!=std::string::npos){
    s.erase(pos,1);
  }
}**/


bool msgfdbReader::checkValidity(const std::string file)
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
  
  if (line.find("SpecIndex") != std::string::npos && line.find("MSGFScore") != std::string::npos) //NOTE there doesn't seem to be any good way to check if the file is from msgfdb
  {
    std::vector<std::string> column_names=split(line,'\t');
    if(!(column_names.size()==14||column_names.size()==15))//Check that the size is corrrect, it should have length 14 or 15
    { 
      std::cerr << "The file " << file << " has the wrong number of columns: "<< column_names.size() << ". Should be 14 or 15." << std::endl;
      exit(-1);
    }
  } 
  else
  {
    ismeta = false;
  }

  std::cerr << "TMP validate"<< std::endl;
  
  
  
  return ismeta;
}


void  msgfdbReader::addFeatureDescriptions(bool doEnzyme,const std::string& aaAlphabet,std::string fn)
{
  
  std::cerr << "TMP featuer decr"<< std::endl;
  
  push_backFeatureDescription("deNovoScore");
  push_backFeatureDescription("MSGFScore");
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
  
  //NOTE NSM?
  //push_backFeatureDescription("lnNumSP");
  
  //Mass difference
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
  
  if (!aaAlphabet.empty()) 
  {
    for (std::string::const_iterator it = aaAlphabet.begin(); it != aaAlphabet.end(); it++)
      push_backFeatureDescription(*it + "-Freq");
  }
}

void msgfdbReader::getMaxMinCharge(const std::string fn){
  
  int n = 0, charge = 0;
  
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  
  getline(sqtIn, line); //First line is column names which is of no intrest for max and min charge

  while (getline(sqtIn, line)) {
    //Get line and look for min/max charge, charge is at pos charge_pos
    std::vector<std::string> psm_vector=split(line,'\t');
    charge=atoi(psm_vector.at(6).c_str()); //TODO Error handling, would be nice if it could handle scientific notation
    if (minCharge > charge) minCharge = charge;
    if (maxCharge < charge) maxCharge = charge;
    n++;
  }
  if (n <= 0) {
    std::cerr << "The file " << fn << " does not contain any records"<< std::endl;
    sqtIn.close();
    exit(-1);
  }

  std::cerr << "TMP search"<< std::endl;
  
  
  sqtIn.close();
  return;
}

void msgfdbReader::readPSM(std::string line,bool isDecoy,std::string fileId,
			   boost::shared_ptr<FragSpectrumScanDatabase> database, std::vector<std::string> column_names){
  
  std::ostringstream id;
  std::vector< std::string > proteinIds;
  std::vector< std::string > proteinIdsDecoys; //This vector is only used when we read a combined file
  
  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::map<char,int> ptmMap = po.ptmScheme;
  
  std::vector<std::string> psm_vector=split(line,'\t');
  if(psm_vector.size()!=column_names.size())
  {
    std::cerr << "One row or more in " << fileId << " has the wrong number of columns: "<< psm_vector.size() << ". Should be " << column_names.size() << std::endl;
    exit(-1);
  }
  
  //Variables related to the MSGFDB file type 
  //TODO error handling
  std::string specFile=psm_vector.at(0);
  double specIndex=atof(psm_vector.at(1).c_str());
  double scan=atof(psm_vector.at(2).c_str());
  std::string fragMethod=psm_vector.at(3);
  double observedMassCharge=atof(psm_vector.at(4).c_str()); //Called precursor mass
  double pmError=atof(psm_vector.at(5).c_str());
  int charge=atoi(psm_vector.at(6).c_str()); //TODO would be nice if it could handle scientific notation
  std::string peptide=psm_vector.at(7);
  std::string proteinID=psm_vector.at(8);
  double deNovoScore=atof(psm_vector.at(9).c_str());
  double MSGFScore=atof(psm_vector.at(10).c_str());
  double specProb=atof(psm_vector.at(11).c_str());
  double pValue=atof(psm_vector.at(12).c_str());
  double FDR=atof(psm_vector.at(13).c_str()); //if -tda 1 option was NOT used in msgfdb this is really EFDR specfied in the msgfdb documentation
  if(column_names.at(13)=="FDR"){//-tda 1 option is was used in msgfdb there is one more feature
    double PepFDR=atof(psm_vector.at(14).c_str());
  }
  
  //Create id
  id.str("");
  id << fileId << '_' << specIndex << '_' << scan;
  std::string psmId=id.str();
  
  //Get rid of unprinatables in proteinID and make list of proteinIDs
  proteinID=getRidOfUnprintables(proteinID);
  if(po.iscombined){
    if( (line.find(po.reversedFeaturePattern) != std::string::npos) )
      proteinIdsDecoys.push_back(proteinID);
    else
      proteinIds.push_back(proteinID);
  }
  else{
    proteinIds.push_back(proteinID);
  }
  
  //Check length of peptide, if its to short it cant contain both flanks and peptide
  assert(peptide.size() >= 5 );
  
  //Calculate peptide mass
  double calculatedMassToCharge=Reader::calculatePepMAss(peptide,charge);

  std::cerr << "TMP3  "<< specIndex << std::endl;
  
  //Push_back the DeNovoScore and msgfscore
  f_seq.push_back(deNovoScore);
  f_seq.push_back(MSGFScore);

  f_seq.push_back( observedMassCharge ); // Observed mass
  f_seq.push_back( peptideLength(peptide)); // Peptide length
  int nxtFeat = 8;
  for (int c = minCharge; c
  <= maxCharge; c++)
    f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge

  if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) {
    f_seq.push_back( Enzyme::isEnzymatic(peptide.at(0),peptide.at(2)) ? 1.0 : 0.0);
    f_seq.push_back(Enzyme::isEnzymatic(peptide.at(peptide.size() - 3),peptide.at(peptide.size() - 1)) ? 1.0 : 0.0);
    std::string peptid2 = peptide.substr(2, peptide.length() - 4);
    f_seq.push_back( (double)Enzyme::countEnzymatic(peptid2) );
  }
  
  //Calculate difference between observed and calculated mass
  double dM =MassHandler::massDiff(observedMassCharge, calculatedMassToCharge,
                                   charge, peptide.substr(2, peptide.size()- 4));
  f_seq.push_back( dM ); // obs - calc mass
  f_seq.push_back( (dM < 0 ? -dM : dM)); // abs only defined for integers on some systems
  
  if (po.calcPTMs) 
  {
    f_seq.push_back(cntPTMs(peptide));
  }
  if (po.pngasef) 
  {
    f_seq.push_back(isPngasef(peptide,isDecoy));
  }
  if (po.calcAAFrequencies) {
    computeAAFrequencies(peptide, f_seq);
  }
  
  std::cerr << "TMP remove mod  "<< std::endl;
  
  //Get the flanks/termini and remove them from the peptide sequence
  std::vector<std::string> tmp_vect=split(peptide,'.');
  
  std::cerr << "TMP4 1 "<< peptide << std::endl;
  
  for(int i=0; i<tmp_vect.size();i++){
    std::cerr << "asdf  "<< tmp_vect.at(i) << std::endl;
  }
  
  percolatorInNs::occurence::flankN_type flankN=tmp_vect.at(0);
  
  std::cerr << "TMP4 2 "<< std::endl;
  
  percolatorInNs::occurence::flankC_type flankC=tmp_vect.at(2); //Fel
  
  std::cerr << "TMP4 3 "<< std::endl;
  
  std::string peptideSequence=tmp_vect.at(1);
  
  std::cerr << "TMP4 4 "<< std::endl;
  
  std::string peptideS = peptideSequence;
  
  std::cerr << "TMP 6  "<< std::endl;
  
  //Remove modifications
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
  
  std::cerr << "TMP 7  "<< std::endl;
  
  //Register the ptms
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
  
  std::cerr << "TMP iscombied  "<< std::endl;
  
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

void msgfdbReader::read(const std::string fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database)
{
  std::string fileId;
  
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  
  //The first row contains the names of the columns
  getline(sqtIn, line);
  std::vector<std::string> column_names=split(line,'\t');
  
  if(column_names.at(0)[0]=='#'){
    column_names.at(0)=column_names.at(0).erase(0,1); //Remove "#" in the first name
  }
  
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) fileId.erase(0, spos + 1);
  spos = fileId.find('.');
  if (spos != std::string::npos) fileId.erase(spos);
  
  //Read file line by line, each line is one psm. Tab delimated
  while (getline(sqtIn, line)) { //TODO Fix max psm thing po.hitsPerSpectrum

    std::cerr << "TMP line iter"<< std::endl;

    readPSM(line,isDecoy,fileId, database, column_names);
  }
}

