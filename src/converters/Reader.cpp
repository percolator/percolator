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
  
  bool isMeta = false;
  
  if(po.readProteins)
  {
    readProteins(po.targetDb,po.decoyDb);
  }
  
  //check files are metafiles or not
  if(!po.iscombined)
  {
    bool isMetaTarget = checkIsMeta(po.targetFN); 
    bool isMetaDecoy = checkIsMeta(po.decoyFN);
    if(isMetaTarget == isMetaDecoy)
      isMeta = isMetaTarget;
    else
    {
      std::cerr << "ERROR : one of the input files is a metafile whereas the other one is not a metaFile. " << std::endl;
    }
  }
  else
    isMeta= checkIsMeta(po.targetFN);
  
  //NOTE getMaxMinCharge does more than get max charge and min charge for some types of converters. tandemReader for instance checks if a certain attribute is present or not.
  
  //if they are metaFiles check for max/min charge in the meta of the target files
  //otherwise check for max/min charge in the target file
  if(isMeta)
  {
    //I iterate over all the files to check validity and get max/min charge but I only get the charge from the 
    // target files, max/min charge should be the same for the decoy files
    std::string line;
    std::ifstream meta(po.targetFN.data(), std::ios::in);
    while (getline(meta, line)) {
      if (line.size() > 0 && line[0] != '#') {
	 line.erase(std::remove(line.begin(),line.end(),' '),line.end());
	 checkValidity(line);
	 getMaxMinCharge(line,false);
      }
     }
    meta.close();
    if(!po.iscombined)
    {
      meta.open(po.decoyFN.data(), std::ios::in);
      while (getline(meta, line)) {
	if (line.size() > 0 && line[0] != '#') {
	  line.erase(std::remove(line.begin(),line.end(),' '),line.end());
	  checkValidity(line);
	  getMaxMinCharge(line,true);
	}
      }
      meta.close();
    }
   }
  else
  {
    checkValidity(po.targetFN);
    getMaxMinCharge(po.targetFN,false);
    if(!po.iscombined){
      checkValidity(po.decoyFN);
      getMaxMinCharge(po.decoyFN,true);
    }
  }
  //once I have max/min charge I can put in the features
  addFeatureDescriptions(Enzyme::getEnzymeType() != Enzyme::NO_ENZYME,po.calcAAFrequencies ? aaAlphabet : "");
  
  if (!po.iscombined) 
  {
    std::cerr << "Reading input from files: " <<  std::endl;
    translateFileToXML(po.targetFN, false /* is_decoy */,0,isMeta);
    translateFileToXML(po.decoyFN, true /* is_decoy */,0,isMeta);
  } 
  else 
  {
    
    std::cerr << "Reading input from a combined (target-decoy) file .." << std::endl;
    translateFileToXML(po.targetFN, false /* is_decoy */,0,isMeta);
  }

  // read retention time if the converter was invoked with -2 option
  if (po.spectrumFN.size() > 0) {
    readRetentionTime(po.spectrumFN);
    databases[0]->initRTime(&scan2rt);
    storeRetentionTime(databases[0]);
  }

  xercesc::XMLPlatformUtils::Terminate();
}

void Reader::print(ofstream &xmlOutputStream)
{
  
  xercesc::XMLPlatformUtils::Initialize ();
  
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
  
  if(po.readProteins)
  {    
    serializer ser;
    if (po.xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);

    std::auto_ptr< ::percolatorInNs::parameters > parameters
    ( new ::percolatorInNs::parameters(po.missed_cleavages,po.peptidelength,po.maxmass,po.minmass,po.maxpeplength,po.iscombined));
    ::percolatorInNs::databases databases(po.targetDb,po.decoyDb,parameters);
    
    ser.next ( PERCOLATOR_IN_NAMESPACE, "databases",databases);
  }
  
  string commandLine = "\n<process_info>\n" +
      string("  <command_line>") + po.call.substr(0,po.call.length()-1)
      + "</command_line>\n" + "</process_info>\n";
      
  if (po.xmlOutputFN == "") 
    cout << commandLine;
  else 
    xmlOutputStream << commandLine;

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

  if(po.readProteins && !proteins.empty())
  { 
    if (po.xmlOutputFN == "") 
      //cout << "\n<proteins>\n";
      cout << "\n";
    else 
      //xmlOutputStream << "\n<proteins>\n";
      xmlOutputStream << "\n";
    
    serializer ser;
    std::vector<Protein*>::const_iterator it;
    if (po.xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    for (it = proteins.begin(); it != proteins.end(); it++) 
    { //NOTE uss this is horrible, I should serialize in a Btree the object protein as the PSMs
      //FIXME thi serialization is creating a gap \o between elements
      std::auto_ptr< ::percolatorInNs::protein> p (new ::percolatorInNs::protein((*it)->name,(*it)->length,
							(*it)->totalMass,(*it)->sequence,(*it)->id,(*it)->isDecoy));
      ser.next(PERCOLATOR_IN_NAMESPACE, "protein", *p);
    }
    
    if (po.xmlOutputFN == "") 
      //cout << "\n</proteins>\n";
      cout << "\n";
    else 
      //xmlOutputStream << "\n</proteins>\n";
      xmlOutputStream << "\n";
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


void Reader::translateFileToXML(const std::string fn, bool isDecoy, unsigned int lineNumber_par, bool isMeta)
  {
      
    if(!isMeta)
    {
      // there must be as many databases as lines in the metafile containing the
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

      read(fn,isDecoy,databases[lineNumber_par]);
    
    } else {
      // we hopefully found a meta file
      std::cerr << "Found a meta file" << std::endl;
      unsigned int lineNumber=0;
      std::string line2;
      std::ifstream meta(fn.data(), std::ios::in);
      while (getline(meta, line2)) {
	if (line2.size() > 0 && line2[0] != '#') {
	  line2.erase(std::remove(line2.begin(),line2.end(),' '),line2.end());
	  translateFileToXML(line2,isDecoy,lineNumber,false);
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
  if (pepsequence.length () >= po.peptidelength) {
    
    for(unsigned i=0; i<pepsequence.length();i++)
    {
      if(isalpha(pepsequence[i])){
        mass += massMap_[pepsequence[i]];
      }
    }
    
    mass = (mass + massMap_['o'] + (charge * massMap_['h']) + 1.00727649);
  }
  else
  {
    std::cerr << "Calculate peptide mass error: The peptide is to short: " << pepsequence << std::endl;
    exit(-1);
  }
  
  return(mass);
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


void Reader::parseDataBase(const char* seqfile, bool isDecoy, bool isCombined, unsigned &proteins_counter)
{

  //NOTE I do not parse the gene names because there is not standard definiton for it
  try
  {
    filebuf fb;
    fb.open (seqfile,ios::in);
    istream buffer(&fb);
    std::cerr << "Reading fasta file : " << seqfile << std::endl;
    if (!buffer.eof()) {
	  wchar_t c;
	  while (!buffer.eof()) {
	    c = buffer.peek();
	    if( c == '>' )
	    {
	      std::string protein_name;
	      std::string protein_seq;
	      read_from_fasta(buffer,protein_name,protein_seq);
	      //std::cerr << " Reading " << protein_name << " " << protein_seq << std::endl;
	      if(isCombined) isDecoy = protein_name.find(po.reversedFeaturePattern,0) != std::string::npos;
	      std::set<std::string> peptides;
	      double totalMass = 0.0;
	      //NOTE here I should check the enzyme and do the according digestion
	      //unsigned num_tryptic = calculateProtLengthElastase(protein_seq,peptides,totalMass);
	      unsigned num_tryptic = calculateProtLengthTrypsin(protein_seq,peptides,totalMass);
	      //unsigned num_tryptic = calculateProtLengthChymotrypsin(protein_seq,peptides,totalMass);
	      //unsigned num_tryptic = calculateProtLengthThermolysin(protein_seq,peptides,totalMass);
	      //unsigned num_tryptic = calculateProtLengthProteinasek(protein_seq,peptides,totalMass);
	      Protein *tmp = new Protein();
	      tmp->id = ++proteins_counter;
	      tmp->name = protein_name;
	      tmp->isDecoy = isDecoy;
	      tmp->sequence = protein_seq;
	      tmp->peptides = peptides;
	      tmp->length = num_tryptic;
	      tmp->totalMass = totalMass;
	      proteins.push_back(tmp);
	   }
	 } 
    }
    else
    {
      std::cerr <<  "Error reading combined database : " << seqfile <<  std::endl;
      exit(-1);
    }
  
    buffer.clear();
    fb.close();
    
  }catch(const std::exception &e)
  {
    e.what();
  }
  
  std::string type = isDecoy ?  "target" : "decoy";
  if(po.iscombined) type = "target and decoy";
  
  if(proteins_counter == 0)
  {
    std::cerr << "Error parsing the database, the database given does not contain " << type << " proteins\n" << std::endl;
    exit(-1);
  }
  else if(VERB > 2)
  {
    std::cerr << "\nRead " << proteins_counter << type << " proteins\n" << std::endl;
  }
  
  return;
}

unsigned int Reader::calculateProtLengthTrypsin(string protsequence, 
						set< string >& peptides, double& totalMass)
{
  size_t length = protsequence.length();
  
  if (length>0 && protsequence[length-1]=='*') {
    --length;
  }
  
  for(size_t start=0; start<length; start++)
  {
    if( start == 0 || ( protsequence[start+1] != 'P' && ( protsequence[start] == 'K' || protsequence[start] == 'R' ) ) ) 
    {
      int numMisCleavages = 0;  
      for(size_t end=start+1;( end<length && numMisCleavages <= po.missed_cleavages );end++)
      {
	//NOTE I am missing the case when a tryptip digested peptide has a K|R and P at the end of the sequence
        if( (protsequence[end] == 'K' || protsequence[end] == 'R') && protsequence[end+1] != 'P' )
	{
	   int begin = start;
	   int finish = end - start + 1;
	   if(start != 0)
	   {  
	     begin++;
	     finish--;
	   }
	   std::string peptide = protsequence.substr(begin,finish);
	   double  mass = calculatePepMAss(peptide);
	   
	   if((mass > po.minmass) && (mass< po.maxmass) && (peptide.size() >= po.peptidelength))
	   {
	     peptides.insert(peptide);
	     totalMass+=mass;
	   }
	    numMisCleavages++;
        }
      } 
    }
  }
  
  unsigned size = peptides.size();
  return size;
}

void Reader::readProteins(string filenameTarget, string fileNamedecoy)
{
    unsigned proteins_counter = 0; 
    //NOTE improve the error testing cases
    if(fileNamedecoy == "" && filenameTarget != "")
    {
      parseDataBase(filenameTarget.c_str(),false,true,proteins_counter);
    }
    else if(filenameTarget != "" && fileNamedecoy != "")
    {
      parseDataBase(filenameTarget.c_str(),false,false,proteins_counter);
      parseDataBase(fileNamedecoy.c_str(),true,false,proteins_counter);
    }
    else
    {  
      std::cerr << "\nError database file/s could not be loaded\n" << std::endl;
      exit(-1);
    }
}

void Reader::read_from_fasta(istream &buffer, std::string &name , std::string &seq) 
{
  std::string comment;
  char b[BUFFER_LEN+1];
  b[BUFFER_LEN] = '\0';
  stringstream buf;
  char c = buffer.peek();

   while (!buffer.eof() and (c == ' ' or c == '\n')) {   
      buffer.ignore(1);
      if (!buffer.eof()) c = buffer.peek();
   }
   
   if (!buffer.eof() and c != '>') {
      stringstream ss;
      ss << "next character is " << c;
      std::cerr << "ERROR parsing fasta " << "Incorrect format " << ss.str() << std::endl;
      exit(-1);
   }
   
  buffer.getline(b, BUFFER_LEN);
  name = string(b+1);
  name.erase(std::remove(name.begin(), name.end(), '\r'), name.end());

  long int pos;
  {
    long int pos1 = name.find(' ');
    long int pos2 = name.find('\t');
    if (pos1 != -1 and pos2 != -1)
      pos = pos1 < pos2 ? pos1 : pos2;
    else if (pos1 != -1)
      pos = pos1;
    else if (pos2 != -1)
      pos = pos2;
    else
      pos = -1;
  }

  if (pos > 0) {
    comment = name.substr(pos+1,name.length());
    name = name.substr(0,pos);
  }
  
  string temp;
  char char_temp;
  
  while (!buffer.eof() and (buffer.peek() != '>')) {
    
    buffer >> temp;
    for (string::iterator iter = temp.begin(); iter != temp.end(); iter++) {

      if ((*iter >= 'A') and (*iter <= 'Z'))
      {
	buf << *iter;
      }
      else
      {
	std::cerr << "ERROR parsing fasta " << "Incorrect fasta sequence character " << *iter << std::endl;
	exit(-1);
      }
    }
  
    c = buffer.peek();

    while (!buffer.eof() and (c == ' ' or c == '\n' or c == '\r')) {
    
      buffer.ignore(1);
      if (!buffer.eof()) c = buffer.peek();
    }
 }
 
  seq = string(buf.str());
  
  return;

}