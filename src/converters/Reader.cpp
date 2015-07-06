#include "Reader.h"
#include <typeinfo>

const std::string Reader::aaAlphabet("ACDEFGHIKLMNPQRSTVWY");
const std::string Reader::ambiguousAA("BZJX");
const std::string Reader::modifiedAA("#@*");
const std::string Reader::additionalAA("UO");
const std::string Reader::freqAA(aaAlphabet + ambiguousAA + additionalAA);

//THIS MAP SHOULD BE CREATED FROM THE UNIMOD PTMS xml document
const std::map<unsigned, double> Reader::ptmMass =
     boost::assign::map_list_of(35, 0.0) (21, 0.0) (28, 0.0) (4, 0.0)
                                (1, 0.0) (214, 0.0) (39, 0.0) (7, 0.0)
                                (730, 0.0) (364, 0.0) (29, 0.0) (27, 0.0);


Reader::Reader(ParseOptions *__po)
:po(__po)
{
  tmpDirs = std::vector<char*>();
  tmpFNs = std::vector<std::string>();
  maxCharge = -1;
  minCharge = 10000;
  initMassMap(po->monoisotopic);
}

Reader::Reader() {
  po = 0;
  tmpDirs = std::vector<char*>();
  tmpFNs = std::vector<std::string>();
  maxCharge = -1;
  minCharge = 10000;
  initMassMap(po->monoisotopic);
}


Reader::~Reader() {
  for (unsigned int i=0; i<tmpDirs.size(); i++) {
    if (boost::filesystem::is_directory(tmpDirs[i])) {
      boost::filesystem::remove_all(tmpDirs[i]);
    }
    free(tmpDirs[i]);
  }
}

void Reader::init() {
  // initializing xercesc objects corresponding to pin element...
  xercesc::XMLPlatformUtils::Initialize ();

  // ... <featureDescriptions>
  std::auto_ptr<percolatorInNs::featureDescriptions> fdes_p (new ::percolatorInNs::featureDescriptions());

  // ... <process_info>
  percolatorInNs::process_info::command_line_type command_line = po->call;
  std::auto_ptr<percolatorInNs::process_info> proc_info (new ::percolatorInNs::process_info(command_line));

  // ... <experiment>
  std::auto_ptr< ::percolatorInNs::experiment > ex_p (new ::percolatorInNs::experiment("mitt enzym", proc_info, fdes_p));

  f_seq = ex_p->featureDescriptions();
  fss = ex_p->fragSpectrumScan();

  bool isMeta = false;

  if (po->readProteins) {
    readProteins(po->targetDb,po->decoyDb);
  }

  // check files exists and if they are metafiles or not
  if (!po->iscombined) {
    std::ifstream targetFileIn(po->targetFN.data(), std::ios::in);
    std::ifstream decoyFileIn(po->decoyFN.data(), std::ios::in);
    if (!targetFileIn) {
      targetFileIn.close();
      decoyFileIn.close();
      ostringstream temp;
      temp << "Error : unable to open or read file " << po->targetFN << std::endl;
      throw MyException(temp.str());
    } else if (!decoyFileIn) {
      targetFileIn.close();
      decoyFileIn.close();
      ostringstream temp;
      temp << "Error : unable to open or read file " << po->decoyFN << std::endl;
      throw MyException(temp.str());
    }
    bool isMetaTarget = checkIsMeta(po->targetFN);
    bool isMetaDecoy = checkIsMeta(po->decoyFN);
    if (isMetaTarget == isMetaDecoy) {
      isMeta = isMetaTarget;
    } else {
      ostringstream temp;
      temp << "Error : one of the input files is a metafile whereas"
		       << " the other one is not a metaFile. " << std::endl;
      throw MyException(temp.str());
    }
  } else {
    isMeta = checkIsMeta(po->targetFN);
  }
  // NOTE getMaxMinCharge does more than get max charge and min charge for some types of converters.
  //      tandemReader for instance checks if a certain attribute is present or not.
  if (isMeta) {
    std::string line;
    std::ifstream meta(po->targetFN.data(), std::ios::in);
    while (getline(meta, line)) {
      if (line.size() > 0 && line[0] != '#') {
        line.erase(std::remove(line.begin(),line.end(),' '),line.end());
        checkValidity(line);
        getMaxMinCharge(line,false);
      }
    }
    meta.close();
    if (!po->iscombined) {
      meta.open(po->decoyFN.data(), std::ios::in);
      while (getline(meta, line)) {
        if (line.size() > 0 && line[0] != '#') {
          line.erase(std::remove(line.begin(),line.end(),' '),line.end());
          checkValidity(line);
          getMaxMinCharge(line,true);
        }
      }
      meta.close();
    }
  } else {
    checkValidity(po->targetFN);
    getMaxMinCharge(po->targetFN,false);
    if (!po->iscombined){
      checkValidity(po->decoyFN);
      getMaxMinCharge(po->decoyFN,true);
    }
  }
  //once I have max/min charge I can put in the features
  addFeatureDescriptions(Enzyme::getEnzymeType() != Enzyme::NO_ENZYME);

  if (!po->iscombined) {
    translateFileToXML(po->targetFN, false /* is_decoy */,0,isMeta);
    translateFileToXML(po->decoyFN, true /* is_decoy */,0,isMeta);
  } else {
    translateFileToXML(po->targetFN, false /* is_decoy */,0,isMeta);
  }

  // read retention time if the converter was invoked with -2 option
  if (po->spectrumFN.size() > 0) {
    readRetentionTime(po->spectrumFN);
    databases[0]->initRTime(&scan2rt);
    storeRetentionTime(databases[0]);
  }

  xercesc::XMLPlatformUtils::Terminate();
}

// print to cout or populate xml/tab file
void Reader::print(ostream &outputStream, bool xmlOutput) {  
  if (xmlOutput) {
    xercesc::XMLPlatformUtils::Initialize ();

    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    string headerStr = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n" +
        string("<experiment xmlns=\"") + PERCOLATOR_IN_NAMESPACE + "\"" +
        " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" +
        " xsi:schemaLocation=\"" + PERCOLATOR_IN_NAMESPACE +
        " https://github.com/percolator/percolator/raw/pin-" + schema_major +
        "-" + schema_minor + "/src/xml/percolator_in.xsd\"> \n";

    outputStream << headerStr;
    if (VERB>2 && po->xmlOutputFN == "")
      cerr <<  "The output will be written to " << po->xmlOutputFN << endl;
    
    string enzymeStr = "\n<enzyme>" + Enzyme::getStringEnzyme() + "</enzyme>\n";

    outputStream << enzymeStr;

    if (po->readProteins) {
      serializer ser;
      ser.start (outputStream);

      ::percolatorInNs::databases databases(po->targetDb,po->decoyDb);
      ser.next ( PERCOLATOR_IN_NAMESPACE, "databases",databases);
    }

    string commandLine = "\n<process_info>\n" +
        string("  <command_line>") + po->call.substr(0,po->call.length()-1)
        + "</command_line>\n" + "</process_info>\n";
    outputStream << commandLine;
    
    if (VERB>2)
       cerr << "\nWriting output:\n";
    
    // print features (names + default values)
    { // brackets guarantee that xerces objects go out of scope before Terminate()
      serializer ser;
      ser.start (outputStream);
      ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",f_seq);
    }
    
    // print fragSpectrumScans
    if (VERB>2)
      std::cerr << "Databases : " << databases.size() << std::endl;

    for (unsigned int i=0; i<databases.size();i++) {
      serializer ser;
      ser.start (outputStream);
      if (VERB>2){
        cerr << "outputting content of " << databases[i]->id
            << " (and correspondent decoy file)\n";
      }
      databases[i]->print(ser);
      databases[i]->terminate();
    }
    
    // print proteins
    if (po->readProteins && !proteins.empty()) {
      outputStream << "\n";

      serializer ser;
      std::vector<Protein*>::const_iterator it;
      ser.start (outputStream);
      
      // NOTE I should serialize in a Btree the object protein as the PSMs
      // FIXME the serialization is creating a gap \o between elements
      for (it = proteins.begin(); it != proteins.end(); it++) { 
        std::auto_ptr< ::percolatorInNs::protein > p (new ::percolatorInNs::protein((*it)->name,(*it)->length,
							  (*it)->totalMass,(*it)->sequence,(*it)->id,(*it)->isDecoy));
        ser.next(PERCOLATOR_IN_NAMESPACE, "protein", *p);
      }
      outputStream << "\n";
    }
    
    // print closing tag
    outputStream << "</experiment>" << std::endl;

    xercesc::XMLPlatformUtils::Terminate();
  } else { // tab delimited input files
    if (VERB>2 && po->xmlOutputFN != "")
      cerr <<  "The output will be written to " << po->xmlOutputFN << endl;
    
    if (VERB>2)
      cerr << "\nWriting output:\n";
    
    // print column headers
    bool hasInitialValues = false;
    outputStream << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass";
    BOOST_FOREACH (const ::percolatorInNs::featureDescription & descr, f_seq.featureDescription()) {
      outputStream << "\t" << descr.name();
      if (descr.initialValue().get() != 0) hasInitialValues = true;
    }
    outputStream << "\tPeptide\tProteins" << std::endl;
    
    // print default values
    if (hasInitialValues) {
      outputStream << "DefaultDirection\t-\t-\t-\t-";
      BOOST_FOREACH (const ::percolatorInNs::featureDescription & descr, f_seq.featureDescription()) {
        outputStream << "\t" << descr.initialValue().get();
      }
      outputStream << std::endl;
    }
    
    // print fragSpectrumScans
    if (VERB > 2)
      std::cerr << "Databases : " << databases.size() << std::endl;

    for (unsigned int i = 0; i < databases.size(); i++) {
      databases[i]->printTab(outputStream);
      databases[i]->terminate();
    }
  }
}

void Reader::translateFileToXML(const std::string &fn, bool isDecoy, 
                                unsigned int lineNumber_par, bool isMeta) {
  if (!isMeta) {
    // there must be as many databases as lines in the metafile containing the
    // files. If this is not the case, add a new one
    if (databases.size() == lineNumber_par) {
	    // initialize database
	    std::auto_ptr<serialize_scheme> database(new serialize_scheme(fn));

	    //NOTE this is actually not needed in case we compile with the boost-serialization scheme
	    //indicate this with a flag and avoid the creating of temp files when using boost-serialization
	    if (database->toString() != "FragSpectrumScanDatabaseBoostdb") {
	      // create temporary directory to store the pointer to the database
	      string tcf = "";
	      char * tcd;
	      string str;

	      //TODO it would be nice to somehow avoid these declararions and therefore avoid the linking to
	      //boost filesystem when we dont use them
	      try {
	        boost::filesystem::path ph = boost::filesystem::unique_path();
	        boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
	        boost::filesystem::path file("converters-tmp.tcb");
	        tcf = std::string((dir / file).string());
	        str =  dir.string();
	        tcd = new char[str.size() + 1];
	        std::copy(str.begin(), str.end(), tcd);
	        tcd[str.size()] = '\0';
	        if (boost::filesystem::is_directory(dir)) {
	          boost::filesystem::remove_all(dir);
	        }

	        boost::filesystem::create_directory(dir);
	      } catch (boost::filesystem::filesystem_error &e) {
	        std::cerr << e.what() << std::endl;
	      }

	      tmpDirs.resize(lineNumber_par+1);
	      tmpDirs[lineNumber_par]=tcd;
	      tmpFNs.resize(lineNumber_par+1);
	      tmpFNs[lineNumber_par]=tcf;
	      database->init(tmpFNs[lineNumber_par]);
	    } else {
	      database->init("");
	    }
	    databases.resize(lineNumber_par+1);
	    databases[lineNumber_par]=database;
	    assert(databases.size()==lineNumber_par+1);
    }
    if (VERB>1) {
    	std::cerr << "Reading " << fn << std::endl;
    }

    read(fn,isDecoy,databases[lineNumber_par]);
  } else {
    // we hopefully found a meta file
    if (VERB>1)
      std::cerr << "Found a meta file: " << fn <<std::endl;
      
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
  
void Reader::push_backFeatureDescription(const char * str, const char *description, double initvalue) {
  percolatorInNs::featureDescriptions::featureDescription_sequence &fd_sequence = f_seq.featureDescription();
  std::auto_ptr< ::percolatorInNs::featureDescription > f_p( new ::percolatorInNs::featureDescription(str));
  //adds initial value and description to the description object
  f_p->initialValue(initvalue);
  f_p->description(description);
  assert(f_p.get());
  fd_sequence.push_back(f_p);
  return;
}

string Reader::getRidOfUnprintables(const string &inpString) {
  string outputs = "";
  for (unsigned int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    if (((int)ch) >= 32 && ((int)ch) <= 128) {
      outputs += ch;
    }
  }
  return outputs;
}

void Reader::computeAAFrequencies(const string& pep,  percolatorInNs::features::feature_sequence & f_seq ) {
  //the peptide has to include the flanks
  assert(checkPeptideFlanks(pep));
  // Overall amino acid composition features
  string::size_type aaSize = freqAA.size();

  std::vector< double > doubleV;
  for (unsigned int m = 0  ; m < aaSize ; m++ )  {
    doubleV.push_back(0.0);
  }
  int len = 0;
  for (string::const_iterator it = pep.begin() + 2; it != pep.end() - 2; it++) {
    string::size_type pos = freqAA.find(*it);
    if (pos != string::npos) doubleV[pos]++;
    len++;
  }
  assert(len>0);
  for (unsigned int m = 0  ; m < aaSize ; m++ )  {
    doubleV[m] /= len;
  }
  std::copy(doubleV.begin(), doubleV.end(), std::back_inserter(f_seq));
}

double Reader::calculatePepMAss(const std::string &pepsequence,double charge) {
  double mass  =  0.0;
  assert(!checkPeptideFlanks(pepsequence));

  for(unsigned i=0; i<pepsequence.length();i++) {
    if (freqAA.find(pepsequence[i]) != string::npos) {
      mass += massMap_[pepsequence[i]];
    } else if(modifiedAA.find(pepsequence[i]) != std::string::npos) {
      unsigned annotation = po->ptmScheme[pepsequence[i]];
      mass += ptmMass.at(annotation);
    } else {
      ostringstream temp;
      temp << "Error: estimating peptide mass, the amino acid "
           << pepsequence[i] << " is not valid." << std::endl;
      throw MyException(temp.str());
    }
  }

  mass = (mass + massMap_['o'] + (charge * massMap_['h']) + 1.00727649);
  return mass;
}

std::string Reader::removePTMs(const string& peptide, std::map<char,int>& ptmMap) {
  std::string peptideSequence = peptide;
  if (checkPeptideFlanks(peptide)) {
    peptideSequence = peptide.substr(2, peptide.size()- 4);
  }
  for (unsigned int ix = 0; ix < peptideSequence.size(); ++ix) {
    if (freqAA.find(peptideSequence[ix]) == string::npos) {
      if (peptideSequence[ix] == '[') {
        unsigned int posEnd = peptideSequence.substr(ix).find_first_of(']');
        if (posEnd == string::npos) {
          ostringstream temp;
	        temp << "Error : Peptide sequence " << peptide << " contains an invalid modification" << endl;
          throw MyException(temp.str());
        } else {
          peptideSequence.erase(ix--, posEnd + 1);
        }
      } else if (ptmMap.count(peptideSequence[ix]) > 0) {
	      peptideSequence.erase(ix--,1);
      } else {
        ostringstream temp;
	      temp << "Error : Peptide sequence " << peptide << " contains modification " 
	      << peptideSequence[ix] << " that is not specified by a \"-p\" argument" << endl;
        throw MyException(temp.str());
      }
    }  
  }
  if (checkPeptideFlanks(peptide)) {
    return peptide.substr(0,1) + std::string(".") + peptideSequence + std::string(".") + peptide.substr(peptide.size() - 1,1);
  } else {
    return peptideSequence;
  }
}

unsigned int Reader::peptideLength(const string& pep) {
  unsigned int len = 0;
  assert(checkPeptideFlanks(pep));
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (freqAA.find(pep.at(pos)) != string::npos) {
      len++;
    } else if (pep.at(pos) == '[') {
      unsigned int posEnd = pep.substr(pos).find_first_of(']');
      if (posEnd == string::npos) {
        ostringstream temp;
        temp << "Error : Peptide sequence " << pep << " contains an invalid modification" << endl;
        throw MyException(temp.str());
      } else {
        pos += posEnd;
      }
    }
  }
  return len;
}

unsigned int Reader::cntPTMs(const string& pep, std::map<char,int>& ptmMap) {
  unsigned int len = 0;
  assert(checkPeptideFlanks(pep));
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (ptmMap.find(pep.at(pos)) != ptmMap.end()) {
      ++len;
    } else if (pep.at(pos) == '[') {
      unsigned int posEnd = pep.substr(pos).find_first_of(']');
      if (posEnd == string::npos) {
        ostringstream temp;
        temp << "Error : Peptide sequence " << pep << " contains an invalid modification" << endl;
        throw MyException(temp.str());
      } else {
        ++len;
        pos += posEnd;
      }
    }
  }
  return len;
}

//NOTE this should return bool and the converts to double
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


void Reader::parseDataBase(const char* seqfile, bool isDecoy, bool isCombined, unsigned &proteins_counter) {
  try {
    filebuf fb;
    fb.open (seqfile,ios::in);
    istream buffer(&fb);
    if (VERB>1)
      std::cerr << "Reading fasta file : " << seqfile << std::endl;
    
    if (!buffer.eof()) {
	    wchar_t c;
	    while (!buffer.eof()) {
	      c = buffer.peek();
	      if (c == '>') {
	        std::string protein_name;
	        std::string protein_seq;
	        read_from_fasta(buffer,protein_name,protein_seq);
	        //std::cerr << " Reading " << protein_name << " " << protein_seq << std::endl;
	        if(isCombined) isDecoy = protein_name.find(po->reversedFeaturePattern,0) != std::string::npos;
	        std::set<std::string> peptides;
	        double totalMass = 0.0;
	        //NOTE here I should check the enzyme and do the according digestion a switch
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
    } else {
      fb.close();
      ostringstream temp;
      temp <<  "Error : reading combined database : " << seqfile <<  std::endl;
      throw MyException(temp.str());
    }

    buffer.clear();
    fb.close();

  } catch(const std::exception &e) {
    throw MyException(std::string(e.what()));
  }

  std::string type = isDecoy ?  "target" : "decoy";
  if (po->iscombined) type = "target and decoy";

  if (proteins_counter == 0) {
    ostringstream temp;
    std::cerr << "Error : parsing the database, the database given does not contain "
    << type << " proteins\n" << std::endl;
    throw MyException(temp.str());
  } else if(VERB > 2) {
    std::cerr << "\nRead " << proteins_counter << type << " proteins\n" << std::endl;
  }

  return;
}

unsigned int Reader::calculateProtLengthTrypsin(const string &protsequence,
						set< string >& peptides, double& totalMass) {
  size_t length = protsequence.length();

  if (length>0 && protsequence[length-1]=='*') {
    --length;
  }

  for (size_t start=0; start<length; start++) {
    if (start == 0 || ( protsequence[start+1] != 'P'
        && ( protsequence[start] == 'K' || protsequence[start] == 'R' ) ) ) {
      unsigned int numMisCleavages = 0;
      for(size_t end=start+1;( end<=length && numMisCleavages <= po->missed_cleavages );end++) {
    	//NOTE I am missing the case when a tryptip digested peptide has a K|R and P at the end of the sequence
        if ( (protsequence[end] == 'K' || protsequence[end] == 'R')
	          && (protsequence[end+1] != 'P' || end == length ) ) {
          int begin = start + 1;
          int finish = end - start;
          if (start == 0) {
            begin--;
            finish++;
          }
          std::string peptide = protsequence.substr(begin,finish);
          if (peptide.size() >= po->peptidelength && peptide.size() <= po->maxpeplength) {
	          double mass = calculatePepMAss(peptide);

	          if ((mass > po->minmass) && (mass < po->maxmass) ) {
		          peptides.insert(peptide);
		          totalMass+=mass;
	          }
	        }
          numMisCleavages++;
        }
      }
    }
  }

  return (unsigned)peptides.size();
}

void Reader::readProteins(const string &filenameTarget,const string &fileNamedecoy) {
  unsigned proteins_counter = 0;
  //NOTE improve the error testing cases
  if (fileNamedecoy == "" && filenameTarget != "") {
    parseDataBase(filenameTarget.c_str(),false,true,proteins_counter);
  } else if(filenameTarget != "" && fileNamedecoy != "") {
    parseDataBase(filenameTarget.c_str(),false,false,proteins_counter);
    parseDataBase(fileNamedecoy.c_str(),true,false,proteins_counter);
  } else {
    ostringstream temp;
    std::cerr << "\nError : database file/s could not be loaded\n" << std::endl;
    throw MyException(temp.str());
  }
}

void Reader::read_from_fasta(istream &buffer, std::string &name , std::string &seq) {
  std::string comment;
  char b[BUFFER_LEN+1];
  b[BUFFER_LEN] = '\0';
  stringstream buf;
  char c = buffer.peek();

  while (!buffer.eof() && (c == ' ' || c == '\n')) {
    buffer.ignore(1);
    if (!buffer.eof()) c = buffer.peek();
  }

  if (!buffer.eof() && c != '>') {
    stringstream ss;
    ss << "next character is " << c;
    ss << "Error : parsing fasta file " << "Incorrect format " << ss.str() << std::endl;
    throw MyException(ss.str());
  }

  buffer.getline(b, BUFFER_LEN);
  name = string(b+1);
  name.erase(std::remove(name.begin(), name.end(), '\r'), name.end());

  long int pos;
  {
    long int pos1 = name.find(' ');
    long int pos2 = name.find('\t');
    if (pos1 != -1 && pos2 != -1)
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

  while (!buffer.eof() && (buffer.peek() != '>')) {
    buffer >> temp;
    for (string::iterator iter = temp.begin(); iter != temp.end(); iter++) {
      if ( ((*iter >= 'A') && (*iter <= 'Z')) || (modifiedAA.find(*iter) != std::string::npos) ) {
	      buf << *iter;
      } else {
      	ostringstream temp;
	      temp << "Error : parsing fasta file " << "Incorrect fasta sequence character " << *iter << std::endl;
	      throw MyException(temp.str());
      }
    }

    c = buffer.peek();

    while (!buffer.eof() && (c == ' ' || c == '\n' || c == '\r')) {
      buffer.ignore(1);
      if (!buffer.eof()) c = buffer.peek();
    }
  }
  seq = string(buf.str());

  return;

}

double Reader::massDiff(double observedMass, double calculatedMass,unsigned int charge) {
  assert(charge > 0);
  double dm = observedMass - calculatedMass;
  if (po->monoisotopic) {
    double isodm = dm - 1;
    for (int isotope = 0; isotope < 5; ++isotope) {
      if (abs(isodm) > abs(dm + isotope)) {
        isodm = dm + isotope;
      }
    }
    dm = isodm / calculatedMass;
    return dm;
  }
  return dm / charge;
}

bool Reader::checkPeptideFlanks(const std::string &pep) {
  assert(pep.size() >= 4 );
  if(pep.at(1) == '.' && pep.at(pep.size() - 2) == '.')
    return true;
  else
    return false;
}

void Reader::initMassMap(bool useAvgMass) {
  if (useAvgMass) { // avg masses
    massMap_['h'] =  1.00794;
    massMap_['o'] = 15.9994;
    massMap_['G'] = 57.05192;
    massMap_['A'] = 71.07880;
    massMap_['S'] = 87.07820;
    massMap_['P'] = 97.11668;
    massMap_['V'] = 99.13256;
    massMap_['T'] = 101.10508;
    massMap_['C'] = 103.13880;
    massMap_['L'] = 113.15944;
    massMap_['I'] = 113.15944;
    massMap_['N'] = 114.10384;
    massMap_['D'] = 115.08860;
    massMap_['Q'] = 128.13072;
    massMap_['K'] = 128.17408;
    massMap_['E'] = 129.11548;
    massMap_['M'] = 131.19256;
    massMap_['H'] = 137.14108;
    massMap_['F'] = 147.17656;
    massMap_['R'] = 156.18748;
    massMap_['Y'] = 163.17596;
    massMap_['W'] = 186.21320;

    massMap_['U'] = 150.0588;
    massMap_['O'] = 237.301;
    massMap_['X'] = massMap_['L'];  /* treat X as L or I for no good reason */
    massMap_['B'] = (massMap_['N'] + massMap_['D']) / 2.0;  /* treat B as average of N and D */
    massMap_['Z'] = (massMap_['Q'] + massMap_['E']) / 2.0;  /* treat Z as average of Q and E */
  } else { /* monoisotopic masses */
    massMap_['h'] =  1.0078250;
    massMap_['o'] = 15.9949146;
    massMap_['A'] = 71.0371136;
    massMap_['C'] = 103.0091854;
    massMap_['D'] = 115.0269428;
    massMap_['E'] = 129.0425928;
    massMap_['F'] = 147.0684136;
    massMap_['G'] = 57.0214636;
    massMap_['H'] = 137.0589116;
    massMap_['I'] = 113.0840636;
    massMap_['K'] = 128.0949626;
    massMap_['L'] = 113.0840636;
    massMap_['M'] = 131.0404854;
    massMap_['N'] = 114.0429272;
    massMap_['P'] = 97.0527636;
    massMap_['Q'] = 128.0585772;
    massMap_['R'] = 156.1011106;
    massMap_['S'] = 87.0320282;
    massMap_['T'] = 101.0476782;
    massMap_['U'] = 149.90419;
    massMap_['V'] = 99.0684136;
    massMap_['W'] = 186.07931;
    massMap_['Y'] = 163.06333;

    massMap_['U'] = 150.9536;
    massMap_['O'] = 237.147;
    massMap_['X'] = massMap_['L'];  /* treat X as L or I for no good reason */
    massMap_['B'] = (massMap_['N'] + massMap_['D']) / 2.0;  /* treat B as average of N and D */
    massMap_['Z'] = (massMap_['Q'] + massMap_['E']) / 2.0;  /* treat Z as average of Q and E */
  }
}

void Reader::readRetentionTime(const std::string &filename) {
  MSReader r;
  Spectrum s;
  r.setFilter(MS2);
  char* cstr = new char[filename.size() + 1];
  strcpy(cstr, filename.c_str());
  // read first spectrum
  r.readFile(cstr, s);
  while (s.getScanNumber() != 0) {
    // check whether an EZ line is available
    if (s.sizeEZ() != 0) {
      // for each EZ line (each psm)
      for(int i = 0; i<s.sizeEZ(); i++) {
        // save experimental mass and retention time
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).mh);
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).pRTime);
      }
    } else if((double)s.getRTime() != 0) { // if no EZ line is available, check for an RTime lines
      scan2rt[s.getScanNumber()].push_back(s.getRTime());
    } else { // if neither EZ nor I lines are available
      delete[] cstr;
      ostringstream temp;
      temp << "Error : The ms2 in input file does not appear to contain retention time "
          << "information. Please run without -2 option." << std::endl;
      throw MyException(temp.str());
    }
    // read next scan
    r.readFile(NULL, s);
  }
  delete[] cstr;
}

void Reader::storeRetentionTime(boost::shared_ptr<FragSpectrumScanDatabase> database) {
  // for each spectra from the ms2 file
  typedef std::map<int, vector<double> > map_t;
  BOOST_FOREACH (map_t::value_type& i, scan2rt) {
    // scan number
    int scanNr = i.first;
    // related retention times
    vector<double>* rTimes = &(i.second);
    if (database->getFSS(scanNr).get()!=0) {
      fragSpectrumScan fss = *(database->getFSS(scanNr));
      fragSpectrumScan::peptideSpectrumMatch_sequence& psmSeq = fss.peptideSpectrumMatch();
      // retention time to be stored
      double storeMe = 0;
      // if rTimes only contains one element
      if (rTimes->size()==1) {
        // take that as retention time
        storeMe = rTimes->at(0);
      } else {
        // else, take retention time of psm that has observed mass closest to
        // theoretical mass (smallest massDiff)
        double massDiff = (std::numeric_limits<double>::max)(); // + infinity
        for (fragSpectrumScan::peptideSpectrumMatch_iterator psmIter_i = psmSeq.begin(); psmIter_i != psmSeq.end(); ++psmIter_i) {
          // skip decoy
          if (psmIter_i->isDecoy() != true) {
            double cm = psmIter_i->calculatedMass();
            double em = psmIter_i->experimentalMass();
            // if a psm with observed mass closer to theoretical mass is found
            if (abs(cm-em) < massDiff) {
              // update massDiff
              massDiff = abs(cm-em);
              // get corresponding retention time
              vector<double>::const_iterator r = rTimes->begin();

              // Loop over alternatives EZ-lines, choose the one with the smallest mass difference
              double altMassDiff = (std::numeric_limits<double>::max)();  // + infinity
              for(; r<rTimes->end(); r=r+2) { // Loops over the EZ-line mh values (rounded to one or two decimals...)
                double rrr = *r;  //mass+h
                double exm = psmIter_i->experimentalMass();  //actually masstocharge
                //FIXME: as rrr is m+h and exm is m/z, this ugly fix loops through many charges
                double rrr_mz;
                for(int charge = 1; charge<7; charge++) {
                  rrr_mz = (rrr + (charge-1)*1.007276466) / charge;
                  if(abs(rrr_mz-exm) < altMassDiff) {
                    altMassDiff = abs(rrr_mz-exm);
                    storeMe = *(r+1);
                  }
                }
                std::cerr << rrr << " " << exm << " " << storeMe << " " << altMassDiff << std::endl;
              }
              r = rTimes->end();
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
