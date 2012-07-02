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


#include "msgfdbReader.h"

//In constructor
//Set file name
//Set constants for the psm_array
//Set featureDescriptions


msgfdbReader::msgfdbReader()
{

}

msgfdbReader::msgfdbReader(const msgfdbReader& other)
{

}

msgfdbReader::~msgfdbReader()
{

}

msgfdbReader& msgfdbReader::operator=(const msgfdbReader& other)
{
return *this;
}

bool msgfdbReader::operator==(const msgfdbReader& other) const
{
///TODO: return ...;
}


//Some functions to handle strings
std::vector<std::string> &msgfdbReader::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
      if(item.find_first_of("abcdefghijklmnopqrstuvwxyz0123456789")!=std::string::npos){
	elems.push_back(item); 
      }
    }
    return elems;
}

std::vector<std::string> msgfdbReader::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

std::string msgfdbReader::remove_endl(std::string s)
{
  int pos=s.find('\n');
  if(pos!=std::string::npos){
    s.erase(pos,1);
  }      
}

void msgfdbReader::readPSM(std::string line,bool isDecoy,
    const ParseOptions& po,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss,
     int minCharge, int maxCharge, std::string fileId,
     FragSpectrumScanDatabase* database){
  
  std::ostringstream id;
  std::vector< std::string > proteinIds;
  std::vector< std::string > proteinIdsDecoys; //This vector is only used when we read a combined file
  
  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::map<char,int> ptmMap = po.ptmScheme; // This map should not be const declared as we use operator[], this is a FIXME for frther releases 
  
  //Variables related to the MSGFDB file type TODO error handling
  vector<std::string> psm_vector=split(line,'\t');
  std::string specFile=psm_vector.at(0);
  int specIndex=atoi(psm_vector.at(1).c_str());
  unsigned int scan=atoi(psm_vector.at(2).c_str());
  std::string fragMethod=psm_vector.at(3);
  double observedMassCharge=atof(psm_vector.at(4).c_str()); //Callled precursor mass
  double pmError=atof(psm_vector.at(5).c_str());
  int charge=atoi(psm_vector.at(6).c_str());
  std::string peptide=psm_vector.at(7);
  std::string proteinID=psm_vector.at(8);
  double deNovoScore=atof(psm_vector.at(9).c_str());
  double MSGFScore=atof(psm_vector.at(10).c_str());
  double specProb=atof(psm_vector.at(11).c_str());
  double pValue=atof(psm_vector.at(12).c_str());
  double EFDR=atof(psm_vector.at(13).c_str());
  
  //Create id
  id.str("");
  id << fileId << '_' << specIndex << '_' << scan;
  
  //If combined check if the psm is a decoy or not NOTE What to use it for?
  if(specFile.find(po.reversedFeaturePattern, 0) != std::string::npos && po.iscombined){
    isDecoy=false;
  }else if(po.iscombined){
    isDecoy=true;
  }
  
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
  
  //Calculate peptide mass NOTE Maybee peptide should be without flanks
  MassHandler::initMassMap();
  double calculatedMassToCharge=MassHandler::calculatePepMAss(peptide,charge);

  
  //Push_back the DeNovoScore and msgfscore
  f_seq.push_back(deNovoScore);
  f_seq.push_back(MSGFScore);

  f_seq.push_back( observedMassCharge ); // Observed mass
  f_seq.push_back( DataSet::peptideLength(peptide)); // Peptide length
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
  
  if (po.calcPTMs) f_seq.push_back(  DataSet::cntPTMs(peptide));
  if (po.pngasef) f_seq.push_back( DataSet::isPngasef(peptide, isDecoy));
  if (po.calcAAFrequencies) {
    computeAAFrequencies(peptide, f_seq);
  }
  
  //Get the flanks/termini and remove them from the peptide sequence
  vector<string> tmp_vect=split(peptide,'.');
  percolatorInNs::occurence::flankN_type flankN=tmp_vect.at(0);
  percolatorInNs::occurence::flankC_type flankC=tmp_vect.at(2);
  std::string peptideSequence=tmp_vect.at(1);
  std::string peptideS = peptideSequence;
  
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
  
  if(!po.iscombined)
  {
    
    percolatorInNs::peptideSpectrumMatch* tmp_psm = new percolatorInNs::peptideSpectrumMatch (
	features_p,  peptide_p,id, isDecoy, observedMassCharge, calculatedMassToCharge, charge);
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
	features_p,  peptide_p,id, true /*is decoy*/, observedMassCharge, calculatedMassToCharge, charge);
      
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
	features_p2,  peptide_p2,id, false /*is decoy*/, observedMassCharge, calculatedMassToCharge, charge);
      
      std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(tmp_psm);

      for ( std::vector< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i ) {
	std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
	psm_p->occurence().push_back(oc_p);
      }
    
      database->savePsm(scan, psm_p);
    }
    
  }
}

void msgfdbReader::searchMaxMinCharge(const std::string fn, int* maxCharge,  int* minCharge){
  
  int n = 0, charge = 0;
  
  const int charge_pos=7; //Would be 6 in an array NOTE Duplicate, should be removed when constructor fixed
  
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream sqtIn;
  sqtIn.open(fn.c_str(), std::ios::in);
  if (!sqtIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  
  getline(sqtIn, line); //First line is column names which is not of intrest for max and min charge
  //NOTE Could do an check here for the file and see if its the right type but the file lacks any good way to identify it.
  while (getline(sqtIn, line)) {
    //Get line and look for min/max charge, charge is at pos charge_pos
    vector<string> psm_vector=split(line,'\t');
    charge=atoi(psm_vector.at(charge_pos-1).c_str()); //TODO Error handling
    if (*minCharge > charge) *minCharge = charge;
    if (*maxCharge < charge) *maxCharge = charge;
    n++;
  }
  if (n <= 0) {
    std::cerr << "The file " << fn << " does not contain any records"<< std::endl;
    sqtIn.close();
    exit(-1);
  }

  sqtIn.close();
  return;
}

void msgfdbReader::read(const std::string fn,
    ::percolatorInNs::featureDescriptions & fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,
      bool isDecoy, const ParseOptions & po, int* maxCharge,  int* minCharge,
      parseType pType, FragSpectrumScanDatabase* database) {

  if ( pType == justSearchMaxMinCharge ) {
    searchMaxMinCharge(fn, maxCharge, minCharge);
    return;
  }
  
  std::cerr << "Reading msgfdb " << fn << std::endl;

  int n = 0, charge = 0;
  std::string fileId;
  
  cerr << "\nread start:\n"; //TMP

  const int charge_pos=7; //Would be 6 in an array NOTE Duplicate
  
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
  
  vector<string> column_names=split(line,'\t');
  
  /*
  if(sizeof(column_names)!=n_columns){ //Check that the size is corrrect
    std::cerr << "The file " << fn << " has the wrong number of columns: "<< n << ". Should be: " << n_columns << std::endl;
    exit(-1);
  }
  */
  
  if(column_names.at(0)[0]=='#'){
    column_names.at(0)=column_names.at(0).erase(0,1); //Remove "#" in the first name
  }
  column_names[sizeof(column_names)-1]=remove_endl(column_names[sizeof(column_names)-1]);
  
  if (fds.featureDescription().size() == 0 ) {
    //TODO add descriptions, need to fixe function first
  }
  
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) fileId.erase(0, spos + 1);
  spos = fileId.find('.');
  if (spos != std::string::npos) fileId.erase(spos);
  
  cerr << "\nPSM:\n"; //TMP
  //Read file line by line, each line is one psm. Tab delimated
  while (getline(sqtIn, line)) { //TODO Fix max psm thing po.hitsPerSpectrum
    //Remove \n from the line
    line=remove_endl(line);

    readPSM(line,isDecoy,po,fsss,*minCharge,*maxCharge,fileId, database);
  }
}

