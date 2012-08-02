#include "msgfdbReader.h"

msgfdbReader::msgfdbReader(ParseOptions *po):Reader(po)
{
  
}

msgfdbReader::~msgfdbReader()
{

}


//A function to split strings
std::vector<std::string> &msgfdbReader::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
      if(item.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_#@*?")!=std::string::npos){
	elems.push_back(item); 
      }
    }
    return elems;
}

std::vector<std::string> msgfdbReader::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}


//Check validity of the file
bool msgfdbReader::checkValidity(const std::string &file)
{
  bool isvalid = true;
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
  //NOTE there doesn't seem to be any good way to check if the file is from msgfdb
  if (line.find("SpecIndex") != std::string::npos && line.find("MSGF") != std::string::npos) 
  {
    //Check that the size is corrrect, it should have length 13,14 or 15 depending on settings used in msgfdb
    std::vector<std::string> column_names = split(line,'\t');
    if(!(column_names.size()==13||column_names.size()==14||column_names.size()==15))
    { 
      std::cerr << "The file " << file << " has the wrong number of columns: "<< column_names.size() << ". Should be 13,14 or 15 "  
		 << "depending on which msgfdb options were used." << std::endl;
      exit(-1);
    }
    else
    {
      isvalid=true;
    }
  }
  else
  {
    isvalid = false;
  }

  return isvalid;
}

//Check if its a meta file
bool msgfdbReader::checkIsMeta(const std::string &file)
{
  bool ismeta;
  std::string line;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  getline(fileIn, line);
  fileIn.close();
  //NOTE there doesn't seem to be any good way to check if the file is from msgfdb
  if (line.find("SpecIndex") != std::string::npos && line.find("MSGF") != std::string::npos) 
  {
    ismeta=false;
  } 
  else
  {
    ismeta = true;
  }

  return ismeta;
}

void  msgfdbReader::addFeatureDescriptions(bool doEnzyme)
{
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
  
  //Mass difference
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

//Get max and min charge. Also makes an map of peptide with an set of associated proteinIDs
void msgfdbReader::getMaxMinCharge(const std::string &fn, bool isDecoy){
  
  int n = 0, charge = 0;
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream msgfdbIn;
  msgfdbIn.open(fn.c_str(), std::ios::in);
  if (!msgfdbIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  
  getline(msgfdbIn, line); //First line is column names which is of no intrest for max and min charge

  while (getline(msgfdbIn, line)) {
    //Look for min/max charge
    std::vector<std::string> psm_vector = split(line,'\t');
    charge = boost::lexical_cast<int>(psm_vector.at(6)); //Charge is column 7
    if (minCharge > charge) minCharge = charge;
    if (maxCharge < charge) maxCharge = charge;
    n++;
    
    //Make a map of all peptides and which proteins they are in
    std::string peptide = boost::lexical_cast<std::string>(psm_vector.at(7));
    std::string proteinID = boost::lexical_cast<std::string>(psm_vector.at(8));
    proteinID = getRidOfUnprintables(proteinID);
    
    peptideDecoyKey peptideDecoyPair;
    
    if(isDecoy) peptideDecoyPair = std::make_pair(peptide,1);
    if(!isDecoy) peptideDecoyPair = std::make_pair(peptide,0);
    
    if(peptideProteinMap.count(peptideDecoyPair)==0) //Check if already present
    {
      set<std::string> tmpSet;
      tmpSet.insert(proteinID);
      peptideProteinMap[peptideDecoyPair] = tmpSet;
      
    }
    else
    {
      peptideProteinMap[peptideDecoyPair].insert(proteinID);
    }
  }
  if (n <= 0) 
  {
    std::cerr << "The file " << fn << " does not contain any records"<< std::endl;
    msgfdbIn.close();
    exit(-1);
  }

  msgfdbIn.close();
  return;
}

//Reads the psm from line, calculates features and then saves it.
void msgfdbReader::readPSM(const std::string &line,bool isDecoy,const std::string &fileId,
			   boost::shared_ptr<FragSpectrumScanDatabase> database, 
			   std::vector<std::string> column_names, counterMapType &idCounterMap){
  
  std::ostringstream id,keyStream,psmIdentStream;
  set< std::string > proteinIds;
  bool tda1=false, ndef=false, psmUsed=false;
  int rank;
  double calculatedMassToCharge, calculatedMass, dM;
  peptideDecoyKey keyCounter, peptideDecoyPair;
  psmIdentPairType psmIdentPair;
  percolatorInNs::occurence::flankN_type flankN;
  percolatorInNs::occurence::flankC_type flankC;
  std::string peptideSequence, peptideS, peptideNoFlank;
  
  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  std::map<char,int> ptmMap = po->ptmScheme;
  
  std::vector<std::string> psm_vector=split(line,'\t');
  if(psm_vector.size()!=column_names.size())
  {
    std::cerr << "One row or more in " << fileId << " has the wrong number of columns: "
    << psm_vector.size() << ". Should be " << column_names.size() << std::endl;
    exit(-1);
  }
  if(column_names.size()>13)
  {
    ndef=true; //This means that msgfdb was used with default value to the -n optioner or with -tda 1
    if(column_names.at(13)=="FDR") //If this is true then the -tda 1 option was used in msgfb
    {
      tda1=true;
    }
  }
  
  //Variables related to the MSGFDB file type 
  //NOTE not all of these are actually used
  std::string specFile = 		boost::lexical_cast<std::string>(psm_vector.at(0));
  double specIndex =	  		boost::lexical_cast<double>(psm_vector.at(1));
  double scan =			boost::lexical_cast<double>(psm_vector.at(2));
  std::string fragMethod =		boost::lexical_cast<std::string>(psm_vector.at(3));
  double observedMassCharge =		boost::lexical_cast<double>(psm_vector.at(4)); //Called precursor mass in the msgfdb file
  double pmError =			boost::lexical_cast<double>(psm_vector.at(5));
  int charge =				boost::lexical_cast<int>(psm_vector.at(6));
  std::string peptideWithFlank =	boost::lexical_cast<std::string>(psm_vector.at(7));
  std::string proteinID =		boost::lexical_cast<std::string>(psm_vector.at(8));
  double deNovoScore =			boost::lexical_cast<double>(psm_vector.at(9));
  double MSGFScore =			boost::lexical_cast<double>(psm_vector.at(10));
  double specProb =			boost::lexical_cast<double>(psm_vector.at(11));
  double pValue =			boost::lexical_cast<double>(psm_vector.at(12));
  
  double EFDR; 	//Used if tda1 option was used
  double PepFDR;//Used if tda1 option was used
  double FDR; 	//Used if tda1 option was not used in msgfb
  
  
  if(tda1)
  {//if -tda 1 option is was used in msgfdb there is one more feature and FDR is sligtly different
    EFDR =				boost::lexical_cast<double>(psm_vector.at(13));
    PepFDR =				boost::lexical_cast<double>(psm_vector.at(14));
  }
  else if(ndef)
  {
    FDR =				boost::lexical_cast<double>(psm_vector.at(13).c_str());
  }
  
  //Get the flanks/termini and remove them from the peptide sequence
  std::vector<std::string> tmp_vect=split(peptideWithFlank,'.');
  try{
    flankN=tmp_vect.at(0);
    flankC=tmp_vect.at(2);
    //Sometimes the flank is odd letters like ? _ replace them with -
    //NOTE the split function might sometimes messing things upp if the flank is something odd and not in the list in split.
    if(flankN.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ")==std::string::npos)
    {
      flankN="-";
    }
    if(flankC.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ")==std::string::npos)
    {
      flankC="-";
    }

    peptideSequence=tmp_vect.at(1);
    peptideS = peptideSequence;
  }
  catch(exception e){
    std::cerr << "There is a problem with the peptide string: " << peptideWithFlank << " SpecIndex: " << specIndex << std::endl;
    exit(-1);
  }
  
  peptideNoFlank=flankN + " ." + peptideWithFlank + "." + flankC;

  //Get rid of unprinatables in proteinID and make list of proteinIDs
  proteinID = getRidOfUnprintables(proteinID);
  
  psmIdentStream.str("");
  psmIdentStream << specIndex << "_" << peptideWithFlank;
  
  //If its a combined filed isDecoy will be false when this function is called so it has to adjusted based on pattern
  //Also the keys for the maps will be always be false if its combined filed
  if(po->iscombined)
  {
    isDecoy = proteinID.find(po->reversedFeaturePattern, 0) != std::string::npos;
    //Key for the map with peptides and proteins
    peptideDecoyPair = make_pair (peptideWithFlank,0);
    //Key for the map with psms and if the psm has been used or not
    psmIdentPair = make_pair (psmIdentStream.str(),0);
  }
  else
  {
    //Key for the map with peptides and proteins
    if(isDecoy) peptideDecoyPair = make_pair (peptideWithFlank,1);
    if(!isDecoy) peptideDecoyPair = make_pair (peptideWithFlank,0);
    
    //Key for the map with psms and if the psm has been used or not
    if(isDecoy) psmIdentPair = make_pair (psmIdentStream.str(),1);
    if(!isDecoy) psmIdentPair = make_pair (psmIdentStream.str(),0);
  }
  //If the if statment is not true the function is done beacause this psm, decoy combination has already been put in to the database
  //NOTE is this really needed?
  if(usedPSMs.insert(psmIdentPair).second) 
  {
  
    proteinIds = peptideProteinMap[peptideDecoyPair]; //Get all proteins associated with this peptide
  
    //Get rank from map aka number of hits from the same spectra so far
    keyStream.str("");
    keyStream << specIndex;
    if(isDecoy) keyCounter = make_pair (keyStream.str(),1);
    if(!isDecoy) keyCounter = make_pair (keyStream.str(),0);
    if(idCounterMap.count(keyCounter)==0)
    {
      idCounterMap[keyCounter] = 1;
    }
    else
    {
      idCounterMap[keyCounter]++;
    }
  
    if(rank <= po->hitsPerSpectrum)
    {
      //Create id
      id.str("");
      id << fileId << '_' << specIndex << '_' << charge << '_' << rank;
      std::string psmId = id.str();
  
      //Remove modifications
      for(unsigned int ix=0;ix<peptideSequence.size();++ix) 
      {
	if (aaAlphabet.find(peptideSequence[ix])==string::npos && 
	    ambiguousAA.find(peptideSequence[ix])==string::npos && 
	    additionalAA.find(peptideSequence[ix])==string::npos)
	{
	  if (ptmMap.count(peptideSequence[ix])==0) {
	    cerr << "Peptide sequence " << peptideWithFlank << " contains modification " 
	    << peptideSequence[ix] << " that is not specified by a \"-p\" argument" << endl;
	    exit(-1);
	  }
	  peptideSequence.erase(ix,1);
	}  
      }
      std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSequence   ) );
  
      //Register the ptms
      for(unsigned int ix=0;ix<peptideS.size();++ix) {
	if (aaAlphabet.find(peptideSequence[ix])==string::npos && 
	    ambiguousAA.find(peptideSequence[ix])==string::npos && 
	    additionalAA.find(peptideSequence[ix])==string::npos)
	{
	  int accession = ptmMap[peptideS[ix]];
	  std::auto_ptr< percolatorInNs::uniMod > um_p (new percolatorInNs::uniMod(accession));
	  std::auto_ptr< percolatorInNs::modificationType >  mod_p( new percolatorInNs::modificationType(um_p,ix));
	  peptide_p->modification().push_back(mod_p);      
	  peptideS.erase(ix,1);      
	}  
      }
  
      //Calculate peptide mass
      calculatedMass = calculatePepMAss(peptideNoFlank,charge); 
      //Mass in the msgfdb file is mass/charge
      //NOTE The mass difference is behaving strange, there is somethign wrong somewhere with the mass.
      calculatedMassToCharge = calculatedMass/charge;
  
      //Push_back the DeNovoScore and msgfscore
      f_seq.push_back(deNovoScore);
      f_seq.push_back(MSGFScore);
      f_seq.push_back( observedMassCharge ); // Observed mass
      f_seq.push_back( peptideLength(peptideWithFlank)); // Peptide length
      int nxtFeat = 8;
      for (int c = minCharge; c <= maxCharge; c++)
	f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge

      if (Enzyme::getEnzymeType() != Enzyme::NO_ENZYME) {
	f_seq.push_back( Enzyme::isEnzymatic(peptideWithFlank.at(0),peptideWithFlank.at(2)) ? 1.0 : 0.0);
	f_seq.push_back(Enzyme::isEnzymatic(peptideWithFlank.at(peptideWithFlank.size() - 3),peptideWithFlank.at(peptideWithFlank.size() - 1)) ? 1.0 : 0.0);
	f_seq.push_back( (double)Enzyme::countEnzymatic(peptideNoFlank) );
      }
  
      //Calculate difference between observed and calculated mass
      dM = massDiff(observedMassCharge*charge, calculatedMass, charge);
      f_seq.push_back( dM ); // obs - calc mass
      f_seq.push_back( (dM < 0 ? -dM : dM)); // abs only defined for integers on some systems 
    
      if (po->calcPTMs) 
      {
	f_seq.push_back(cntPTMs(peptideWithFlank));
      }
      if (po->pngasef) 
      {
	f_seq.push_back(isPngasef(peptideWithFlank,isDecoy));
      }
      if (po->calcAAFrequencies) 
      {
	computeAAFrequencies(peptideWithFlank, f_seq);
      }
  
      std::auto_ptr< percolatorInNs::peptideSpectrumMatch >  psm_p(new percolatorInNs::peptideSpectrumMatch (
	features_p,  peptide_p,psmId, isDecoy, observedMassCharge, calculatedMassToCharge, charge));

      for ( set< std::string >::const_iterator i = proteinIds.begin(); i != proteinIds.end(); ++i )
      {
	std::auto_ptr< percolatorInNs::occurence >  oc_p( new percolatorInNs::occurence (*i,flankN, flankC)  );
	psm_p->occurence().push_back(oc_p);
      }
    
      database->savePsm(specIndex, psm_p);

    }//End of if(rank<po.hitsPerSpectrum)
  }//End of if psmUsed
}

void msgfdbReader::read(const std::string &fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database)
{
  std::string fileId;
  std::string line, tmp, prot;
  std::istringstream lineParse;
  std::ifstream msgfdbIn;
  msgfdbIn.open(fn.c_str(), std::ios::in);
  if (!msgfdbIn) {
    std::cerr << "Could not open file " << fn << std::endl;
    exit(-1);
  }
  //The first row contains the names of the columns
  getline(msgfdbIn, line);
  std::vector<std::string> column_names=split(line,'\t');
  
  if(column_names.at(0)[0]=='#'){
    column_names.at(0)=column_names.at(0).erase(0,1); //Remove "#" in the first name
  }
  
  fileId = fn;
  size_t spos = fileId.rfind('/');
  if (spos != std::string::npos) fileId.erase(0, spos + 1);
  spos = fileId.find('.');
  if (spos != std::string::npos) fileId.erase(spos);
  
  if(!idCounterMap.empty())
  {
    idCounterMap.clear();
  }
  //Read file line by line, each line is one psm. Tab delimated
  while (getline(msgfdbIn, line)) {
    readPSM(line,isDecoy,fileId, database, column_names, idCounterMap);
  }
  //Clear map
  idCounterMap.clear();
}

