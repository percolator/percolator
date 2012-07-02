#include "Reader.h"
#include <typeinfo>

const std::string Reader::aaAlphabet("ACDEFGHIKLMNPQRSTVWY");
const std::string Reader::ambiguousAA("BZJX");
const std::string Reader::modifiedAA("#@*");

Reader::Reader(ParseOptions __po)
:po(__po)
{
  tmpDirs = std::vector<char*>();
  tmpFNs = std::vector<std::string>();
  maxCharge = -1;
  minCharge = 10000;
  initMassMap(po.monoisotopic);
}

Reader::~Reader()
{
  //NOTE do this with boost?
  //NOTE remove files as well
  for(int i=0; i<tmpDirs.size(); i++)
    rmdir(tmpDirs[i]);
}


void Reader::init()
{
  // initializing xercesc objects corresponding to pin element...
  xercesc::XMLPlatformUtils::Initialize ();
  
  // ... <featureDescriptions>
  std::auto_ptr<percolatorInNs::featureDescriptions> fdes_p (new ::percolatorInNs::featureDescriptions());

  // ... <process_info>
  percolatorInNs::process_info::command_line_type command_line = po.call;
  std::auto_ptr<percolatorInNs::process_info> proc_info (new ::percolatorInNs::process_info(command_line));

  // ... <experiment>
  std::auto_ptr< ::percolatorInNs::experiment > ex_p (new ::percolatorInNs::experiment("mitt enzym", proc_info, fdes_p));


  f_seq = ex_p->featureDescriptions();
  fss = ex_p->fragSpectrumScan();
  
  if (!po.iscombined) 
  {
    std::cerr << "Reading input from files: " <<  std::endl;
    translateFileToXML(po.targetFN, false /* is_decoy */,0);
    translateFileToXML(po.decoyFN, true /* is_decoy */,0);
  } 
  else 
  {
    
    std::cerr << "Reading input from a combined (target-decoy) file .." << std::endl;
    translateFileToXML(po.targetFN, false /* is_decoy */,0);
  }

  // read retention time if sqt2pin was invoked with -2 option
  if (po.spectrumFN.size() > 0) {
    readRetentionTime(po.spectrumFN);
    databases[0]->initRTime(&scan2rt);
    storeRetentionTime(databases[0]);
  }

  xercesc::XMLPlatformUtils::Terminate();

}

void Reader::print(ofstream &xmlOutputStream)
{
  string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
  string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
  string headerStr = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n" +
      string("<experiment xmlns=\"") + PERCOLATOR_IN_NAMESPACE + "\"" +
      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" +
      " xsi:schemaLocation=\"" + PERCOLATOR_IN_NAMESPACE +
      " https://github.com/percolator/percolator/raw/pin-" + schema_major +
      "-" + schema_minor + "/src/xml/percolator_in.xsd\"> \n";
      
  if (po.xmlOutputFN == "") 
    cout << headerStr;
  else 
  {
    xmlOutputStream << headerStr;
    cerr <<  "The output will be written to " << po.xmlOutputFN << endl;
  }

  string enzymeStr = "\n<enzyme>" + Enzyme::getStringEnzyme() + "</enzyme>\n";
  
  if (po.xmlOutputFN == "") 
    cout << enzymeStr;
  else 
    xmlOutputStream << enzymeStr;

  string commandLine = "\n<process_info>\n" +
      string("  <command_line>") + po.call.substr(0,po.call.length()-1)
      + "</command_line>\n" + "</process_info>\n";
      
  if (po.xmlOutputFN == "") 
    cout << commandLine;
  else 
    xmlOutputStream << commandLine;

  xercesc::XMLPlatformUtils::Initialize ();

  cerr << "\nWriting output:\n";
  // print to cout (or populate xml file)
  // print features
  {
    serializer ser;
    if (po.xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",f_seq);
  }

  // print fragSpecturmScans
  std::cerr << "Databases : " << databases.size() << std::endl;
  for(int i=0; i<databases.size();i++) {
    serializer ser;
    if (po.xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    if(VERB>1){
      cerr << "outputting content of " << databases[i]->id
          << " (and correspondent decoy file)\n";
    }
    databases[i]->print(ser);
    databases[i]->terminte();
  }

  // print closing tag
  if (po.xmlOutputFN == "") 
    std::cout << "</experiment>" << std::endl;
  else 
  {
    xmlOutputStream << "</experiment>" << std::endl;
    xmlOutputStream.close();
  }

  xercesc::XMLPlatformUtils::Terminate();

}


void Reader::translateFileToXML(const std::string fn, bool isDecoy, unsigned int lineNumber_par)
  {
  
    //TODO check its is a metafile or not and if it is valid formed, mmm I should split this function in two to make it clearer
    if(checkValidity(fn))
    {
      // there must be as many databases as lines in the metafile containing sqt
      // files. If this is not the case, add a new one
      if(databases.size()==lineNumber_par)
      {
	// initialize databese     
	std::auto_ptr<serialize_scheme> database(new serialize_scheme(fn));
      
	//NOTE this is actually not needed in case we compile with the boost-serialization scheme
	//indicate this with a flag and avoid the creating of temp files when using boost-serialization
	if(database->toString() != "FragSpectrumScanDatabaseBoostdb")
	{
	  // create temporary directory to store the pointer to the database
	  string tcf = "";
	  char * tcd;
	  string str;
      
	  //TODO it would be nice to somehow avoid these declararions and therefore avoid the linking to
	  //filesystem when we dont use them
	  try
	  {
	    boost::filesystem::path ph = boost::filesystem::unique_path();
	    boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
	    boost::filesystem::path file("converters-tmp.tcb");
	    tcf = std::string((dir / file).string()); 
	    str =  dir.string();
	    tcd = new char[str.size() + 1];
	    std::copy(str.begin(), str.end(), tcd);
	    tcd[str.size()] = '\0';
	    if(boost::filesystem::is_directory(dir))
	    {
	      boost::filesystem::remove_all(dir);
	    }
	
	    boost::filesystem::create_directory(dir);
	  } 
	  catch (boost::filesystem::filesystem_error &e)
	  {
	    std::cerr << e.what() << std::endl;
	  }
	
	  tmpDirs.resize(lineNumber_par+1);
	  tmpDirs[lineNumber_par]=tcd;
	  tmpFNs.resize(lineNumber_par+1);
	  tmpFNs[lineNumber_par]=tcf;
	  database->init(tmpFNs[lineNumber_par]);
	}
	else
	{
	  database->init("");
	}
      
	databases.resize(lineNumber_par+1);
	databases[lineNumber_par]=database;
	assert(databases.size()==lineNumber_par+1);
      }
      if (VERB>1){
	std::cerr << "Reading " << fn << std::endl;
      }
    
      getMaxMinCharge(fn);
      if ( f_seq.featureDescription().size() == 0 ) {
	addFeatureDescriptions(Enzyme::getEnzymeType() != Enzyme::NO_ENZYME,po.calcAAFrequencies ? aaAlphabet : "",po.targetFN);
	//NOTE feature descriptions are the same for different files
      }
      read(fn,isDecoy,databases[lineNumber_par]);
    
    } else {
      // we hopefully found a meta file
      std::cerr << "Found a meta file" << std::endl;
      unsigned int lineNumber=0;
      std::string line2;
      std::ifstream meta(fn.data(), std::ios::in);
      while (getline(meta, line2)) {
	if (line2.size() > 0 && line2[0] != '#') {
	  //NOTE remove the whitespaces
	  line2.erase(std::remove(line2.begin(),line2.end(),' '),line2.end());
	  translateFileToXML(line2,isDecoy,lineNumber);
	  lineNumber++;
	}
      }
      meta.close();
    }
  }
void Reader::push_backFeatureDescription(const char * str) {
  
  percolatorInNs::featureDescriptions::featureDescription_sequence  &fd_sequence =  f_seq.featureDescription();
  std::auto_ptr< ::percolatorInNs::featureDescription > f_p( new ::percolatorInNs::featureDescription(str));
  assert(f_p.get());
  fd_sequence.push_back(f_p);
  return;
}


/**
 * remove non ASCII characters from a string
 */
string Reader::getRidOfUnprintables(string inpString) {
  string outputs = "";
  for (int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    if (((int)ch) >= 32 && ((int)ch) <= 128) {
      outputs += ch;
    }
  }
  return outputs;
}

void Reader::computeAAFrequencies(const string& pep,  percolatorInNs::features::feature_sequence & f_seq ) {
  // Overall amino acid composition features
  assert(pep.size() >= 5);
  string::size_type aaSize = aaAlphabet.size();

  std::vector< double > doubleV;
  for ( int m = 0  ; m < aaSize ; m++ )  {
    doubleV.push_back(0.0);
  }
  int len = 0;
  for (string::const_iterator it = pep.begin() + 2; it != pep.end() - 2; it++) {
    string::size_type pos = aaAlphabet.find(*it);
    if (pos != string::npos) doubleV[pos]++;
    len++;
  }
  assert(len>0);
  for ( int m = 0  ; m < aaSize ; m++ )  {
    doubleV[m] /= len;
  }
  std::copy(doubleV.begin(), doubleV.end(), std::back_inserter(f_seq));
}

double Reader::calculatePepMAss(const std::string &pepsequence,double charge)
{
  double mass  =  0.0;
  if (pepsequence.length () > po.peptidelength) {
    
    for(unsigned i=0; i<pepsequence.length();i++)
    {
      if(isalpha(pepsequence[i])){
        mass += massMap_[pepsequence[i]];
      }
    }
    
    mass = (mass + massMap_['o'] + (charge * massMap_['h'])); 
  }
  return mass; 
}

void Reader::initMassMap(bool useAvgMass)
{
  if (useAvgMass) /*avg masses*/
    {
      massMap_['h']=  1.00794;  
      massMap_['o']= 15.9994;   
      massMap_['G']= 57.05192;
      massMap_['A']= 71.07880;
      massMap_['S']= 87.07820;
      massMap_['P']= 97.11668;
      massMap_['V']= 99.13256;
      massMap_['T']=101.10508;
      massMap_['C']=103.13880;
      massMap_['L']=113.15944;
      massMap_['I']=113.15944;
      massMap_['X']=113.15944;
      massMap_['N']=114.10384;
      massMap_['O']=114.14720;
      massMap_['B']=114.59622;
      massMap_['D']=115.08860;
      massMap_['Q']=128.13072;
      massMap_['K']=128.17408;
      massMap_['Z']=128.62310;
      massMap_['E']=129.11548;
      massMap_['M']=131.19256;
      massMap_['H']=137.14108;
      massMap_['F']=147.17656;
      massMap_['R']=156.18748;
      massMap_['Y']=163.17596;
      massMap_['W']=186.21320;
    }
  else /* monoisotopic masses */
    {
      massMap_['h']=  1.0078250;
      massMap_['o']= 15.9949146;
      massMap_['A']= 71.0371136;
      massMap_['C']=103.0091854;
      massMap_['D']=115.0269428;
      massMap_['E']=129.0425928;
      massMap_['F']=147.0684136;
      massMap_['G']= 57.0214636;
      massMap_['H']=137.0589116;
      massMap_['I']=113.0840636;
      massMap_['K']=128.0949626;
      massMap_['L']=113.0840636;
      massMap_['M']=131.0404854;
      massMap_['N']=114.0429272;
      massMap_['P']= 97.0527636;
      massMap_['Q']=128.0585772;
      massMap_['R']=156.1011106;
      massMap_['S']= 87.0320282;
      massMap_['T']=101.0476782;
      massMap_['U']=149.90419;
      massMap_['V']= 99.0684136;
      massMap_['W']=186.07931;
      massMap_['Y']=163.06333;
    }
}

unsigned int Reader::peptideLength(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (aaAlphabet.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

unsigned int Reader::cntPTMs(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (modifiedAA.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

double Reader::isPngasef(const string& peptide, bool isDecoy ) {
  size_t next_pos = 0, pos;
  while ((pos = peptide.find("N*", next_pos)) != string::npos) {
    next_pos = pos + 1;
    if (! isDecoy) {
      pos += 3;
      if (peptide[pos] == '#') {
        pos += 1;
      }
    } else {
      pos -= 2;
      if (peptide[pos] == '#') {
        pos -= 1;
      }
    }
    if (peptide[pos] == 'T' || peptide[pos] == 'S') {
      return 1.0;
    }
  }
  return 0.0;
}