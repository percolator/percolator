// Erik Sjolund modified this file to handle the percolator xml format
// The original file comes from the streaming example in xsd-3.3.0 ( http://codesynthesis.com/download/xsd/3.3/ )

// original file      : examples/cxx/tree/streaming/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>

#include "mzidentml2pin.h"
#include <boost/lexical_cast.hpp>


typedef map<std::string, mzIdentML_ns::SequenceCollectionType::Peptide_type *> peptideMapType;

typedef map<std::string, int> scanNumberMapType;

std::string call;
enum enzyme_type { enzyme_type_arg_no_enzyme, enzyme_type_arg_elastase, enzyme_type_arg_chymotrypsin, enzyme_type_arg_trypsin};
struct input_options
{
  std::vector<std::string> decoy_file_arg;
  std::vector<std::string> target_file_arg;
  std::string tmp_file_for_indermediate_results_arg;
  bool target_file_given;
  bool decoy_file_given;
  enzyme_type enzyme_type_arg; 
  bool ptm_flag;
  bool pngasef_flag;
  bool aa_freq_flag;

};

static const XMLCh sequenceCollectionStr[] = { chLatin_S, chLatin_e, chLatin_q, chLatin_u, chLatin_e,chLatin_n, chLatin_c, chLatin_e, chLatin_C, chLatin_o, chLatin_l, chLatin_l, chLatin_e, chLatin_c, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };

/* A sketchy overview of the conversion. ( Xpath is used in the explanation )

   We parse the input file(s) two times. The first time ( in function getMinAndMaxCharge() ) is just for finding out the minimum and maximum chargeState. The second time is for all the rest.

   First a hash is created with 

   /mzIdentML/SequenceCollection/Peptide/@id 

   as key and the subtree 

   /mzIdentML/SequenceCollection/Peptide 

   as value.

   Then each 

   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult  

   will be read into memory and translated into a 

   /experiment/fragSpectrumScan

   The first /experiment/fragSpectrumScan/@scan_number will be set to 0, and we increment the @scan_number with +1 for
   the following /experiment/fragSpectrumScan

   We also keep a std::map in memory that maps each 
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/@id 
   to the corresponding /experiment/fragSpectrumScan/@scan_number 

   If we find a /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/@id that is already a key in the std::map
   we don't create a new /experiment/fragSpectrumScan/ but instead merge it into the already created /experiment/fragSpectrumScan/

   The memory consumption of holding all /experiment/fragSpectrumScan/ in memory can be too high so we first store them in a Tokyo Cabinet B+tree database, where we use the @scan_number as key.

   We create our feature descriptions
   /experiment/featureDescriptions/featureDescription 
   by looking at the very first 
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList[0]/SpectrumIdentificationResult[0]/SpectrumIdentificationItem[0]
   and use the 
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList[0]/SpectrumIdentificationResult[0]/SpectrumIdentificationItem[0]/cvParam[ value ]/@name
   as /experiment/featureDescriptions/featureDescription

   Note:
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/SpectrumIdentificationItem/cvParam/@value
   is optional

   so we we restrict with "cvParam[ value ]"
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList[0]/SpectrumIdentificationResult[0]/SpectrumIdentificationItem[0]/cvParam[ value ]/@name

   in other words, we just use cvParam where the attribute "value" is present. In c++ this is "if ( cv.value().present() ) {"

   Each
   /mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/SpectrumIdentificationItem
   translates into a 
   /experiment/fragSpectrumScan/peptideSpectrumMatch

   ---------------------

   The next() method of parser p, needs an explanation:

   It returns:
   first time just the root element
   then the next subtree child of the root element or the next SpectrumIdentificationResult sub tree.
*/

struct MinMaxStruct {
  int min;
  int max;
};

void getMinAndMaxCharge(const char * filename, std::vector< MinMaxStruct > & vec) {
  // This function appends a minMax to the vector
  MinMaxStruct minMax;
  bool foundFirstChargeState = false;
  ifstream ifs;
  ifs.exceptions (ifstream::badbit | ifstream::failbit);
  try
  {
    ifs.open (filename);

    parser p;
    string schemaDefinition = PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, filename, true, schemaDefinition,
      schema_major, schema_minor));
    for (doc = p.next (); doc.get () != 0 && !XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }
    ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());

    int scanNumber = 0;
    for (; doc.get () != 0 && XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan > fss_p( new ::percolatorInNs::fragSpectrumScan( scanNumber, experimentalMassToCharge )); 
      scanNumber++;
      BOOST_FOREACH( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem() )  {
	if ( ! foundFirstChargeState ) {
	  minMax.min = item.chargeState();
	  minMax.max = item.chargeState();
	  foundFirstChargeState = true;
	}
	minMax.min = std::min(item.chargeState(),minMax.min);
	minMax.max = std::max(item.chargeState(),minMax.max);
      }
    }
  }catch (ifstream::failure e) {
    cerr << "Exception opening/reading file :" << filename <<endl;
  }
  catch (const xercesc::DOMException& e)
  {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;  
    XMLString::release(&tmpStr);

  }
  catch (const xml_schema::exception& e)
  {
    cerr << e << endl;
  }
  catch(std::exception e){
    cerr << e.what() <<endl;
  }
  
  assert( foundFirstChargeState );
  vec.push_back( minMax );
  return;
}




void createPSM( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, peptideMapType & peptideMap, int minCharge, int maxCharge, ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge, const input_options & args_info, bool isDecoy, percolatorInNs::featureDescriptions & fdesFirstFile,  ::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_sequence & psm_sequence ) {

  // It is strange but mzIdentML has "experimentalMassToCharge" on the PSM-level in the XML tree. This leads to a lot of redundant information.
  // Let us check that our assumption ( experimentalMassToCharge is constant in the subtree under SpectrumIdentificationResult ) is really valid with an assert()

  assert( item.experimentalMassToCharge() == experimentalMassToCharge );
  assert( item.Peptide_ref().present() );
  assert( peptideMap.find( item.Peptide_ref().get() ) != peptideMap.end() );
  std::string peptideSeq =  peptideMap[item.Peptide_ref().get()]->peptideSequence();
  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSeq ) );
  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  double rank = item.rank();


  if ( ! item.calculatedMassToCharge().present() ) { std::cerr << "error: calculatedMassToCharge attribute is needed for percolator" << std::endl; exit(EXIT_FAILURE); }
  /*
    if (xcorr > 0) {
    f_seq[1] = (xcorr - lastXcorr) / xcorr;
    f_seq[2] = (xcorr - otherXcorr) / xcorr;
    }
    if (!isfinite(f_seq[2])) std::cerr << in;

    f_seq.push_back( log(max(1.0, rSp))); // rank by Sp
  */
  f_seq.push_back( 0.0 ); // delt5Cn (leave until last M line)
  // f_seq.push_back( 0.0 ); // deltCn (leave until next M line)   ................. There is a <cvParam name="sequest:deltacn"
  // f_seq.push_back( xcorr ); There is a  <cvParam name="sequest:xcorr"
  // f_seq.push_back( sp ); There is a <cvParam name="sequest:PeptideRankSp"
  // ???
  //  f_seq.push_back( matched / expected ); // Fraction matched/expected ions    "sequest:matched ions"    "sequest:total ions" 
  //  f_seq.push_back( mass ); // Observed mass
  f_seq.push_back( experimentalMassToCharge * item.chargeState() ); // Observed mass

  // The Sequest mzIdentML format does not have any flankN or flankC so we use "-".
  std::string peptideSeqWithFlanks = std::string("-.") + peptideSeq + std::string(".-");

  f_seq.push_back( DataSet::peptideLength(peptideSeqWithFlanks)); // Peptide length
  int charge = item.chargeState();

  double dM =
    MassHandler::massDiff( item.experimentalMassToCharge(),
			   item.calculatedMassToCharge().get(),
			   charge,
			   peptideSeq );

  for (int c = minCharge; c
	 <= maxCharge; c++) {
    f_seq.push_back( charge == c ? 1.0 : 0.0); // Charge
  }
 
  assert(peptideSeq.size() >= 1 );

  if ( args_info.enzyme_type_arg != enzyme_type_arg_no_enzyme ) {
    f_seq.push_back( Enzyme::isEnzymatic(peptideSeqWithFlanks.at(0),peptideSeqWithFlanks.at(2)) ? 1.0
		     : 0.0);
    f_seq.push_back( 
		    Enzyme::isEnzymatic(peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 3),
					peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 1))
		    ? 1.0
		    : 0.0);
    f_seq.push_back( (double)Enzyme::countEnzymatic(peptideSeq) );
  }
  /*
    f_seq.push_back( log(max(1.0, nSM))); Lukas Kaell told me nSM is not used with Sequest mzIdentML
  */
  f_seq.push_back( dM ); // obs - calc mass
  f_seq.push_back( (dM < 0 ? -dM : dM)); // abs only defined for integers on some systems
  if (args_info.ptm_flag ) { f_seq.push_back(  DataSet::cntPTMs(peptideSeqWithFlanks)); }
  if (args_info.pngasef_flag ) { f_seq.push_back( DataSet::isPngasef(peptideSeqWithFlanks, isDecoy)); }
  //      if (hitsPerSpectrum>1)
  //        feat[nxtFeat++]=(ms==0?1.0:0.0);

  if (args_info.aa_freq_flag ) {
    	  	SqtReader::computeAAFrequencies(peptideSeqWithFlanks, f_seq);
  }
   
  BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_cvParamType & cv, item.cvParam() )  {
    if ( cv.value().present() ) {
      // SpectrumIdentificationItem/cvParam/@value has the datatype string, even though the values seem to be float or double. 
      // percolator_in.xsd uses the datatype double for the features/feature, so we need to convert the string.
      // Using feature_traits for the conversion from std::string to double seems to be the right way to go. Another option would have been to use "strtod()".
      percolatorInNs::features::feature_type fe ( percolatorInNs::features::feature_traits::create(cv.value().get().c_str(),0,0,0) );
      f_seq.push_back( fe);
    }
  }
  BOOST_FOREACH(  const ::mzIdentML_ns::FuGE_Common_Ontology_userParamType & param, item.userParam() )  {
    if ( param.value().present() ) {
      // SpectrumIdentificationItem/userParam/@value has the datatype string, even though the values seem to be float or double or int. 
      // percolator_in.xsd uses the datatype double for the features/feature, so we need to convert the string.
      // Using feature_traits for the conversion from std::string to double seems to be the right way to go. Another option would have been to use "strtod()".
      percolatorInNs::features::feature_type fe ( percolatorInNs::features::feature_traits::create(param.value().get().c_str(),0,0,0) );
      f_seq.push_back( fe);
    }
  }
  // maybe  f_seq.size() normally is greater than  fdesFirstFile.featureDescription().size() ? fix this
  assert( fdesFirstFile.featureDescription().size() == f_seq.size() );
  //is optional for mzIdentML1.0.0 but is compulsory for percolator_in ..... Is this ok?
  assert( item.calculatedMassToCharge().present() );
  ::percolatorInNs::peptideSpectrumMatch* tmp_psm =
      new ::percolatorInNs::peptideSpectrumMatch( features_p, peptide_p, item.id(), isDecoy, item.experimentalMassToCharge(), item.calculatedMassToCharge().get(), item.chargeState());
  std::auto_ptr< ::percolatorInNs::peptideSpectrumMatch > psm_p(tmp_psm);
  psm_sequence.push_back(psm_p);
  return;
}


void addFeatureNameWithEmptyDescription( percolatorInNs::featureDescriptions::featureDescription_sequence & fd_sequence, std::string featureName ) {
  std::auto_ptr< ::percolatorInNs::featureDescription > f_p( new ::percolatorInNs::featureDescription( featureName )); 
  fd_sequence.push_back(f_p);
  return;
}

int loadFromTargetOrDecoyFile( const char * fileName, const input_options & args_info, int minCharge, int maxCharge, bool isDecoy,  percolatorInNs::featureDescriptions & fdesFirstFile, FragSpectrumScanDatabase & database,  int * scanNumber, scanNumberMapType & scanNumberMap  ) {
  namespace xml = xsd::cxx::xml;
  int ret=0;
  
  try
  {
    percolatorInNs::featureDescriptions fdesCurrentFile;
    ifstream ifs;
    ifs.exceptions (ifstream::badbit | ifstream::failbit);
    ifs.open (fileName);
    parser p;
    string schemaDefinition = PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, fileName, true, schemaDefinition,
        schema_major, schema_minor));
    while (doc.get () != 0 && ! XMLString::equals( sequenceCollectionStr, doc->getDocumentElement ()->getTagName())) {
      doc = p.next ();
      // Let's skip some sub trees that we are not interested, e.g. AuditCollection
    }
    assert(doc.get());
    mzIdentML_ns::SequenceCollectionType sequenceCollection(*doc->getDocumentElement ());
    // instead of std::map it should actually be a unordered_map. Let us wait until c++0x finalize.

    peptideMapType peptideMap;

    BOOST_FOREACH( const mzIdentML_ns::SequenceCollectionType::Peptide_type &peptide, sequenceCollection.Peptide() )  {
      assert( peptideMap.find( peptide.id() ) == peptideMap.end() ); // The peptide refs should be unique.
      mzIdentML_ns::SequenceCollectionType::Peptide_type *pept = new mzIdentML_ns::SequenceCollectionType::Peptide_type( peptide);
      assert(pept);
      peptideMap.insert( std::make_pair(peptide.id(), pept ) ) ;
    }
    for (doc = p.next (); doc.get () != 0 && !XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }
    ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
    //      percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence =  ex_p->featureDescriptions().featureDescription();
    percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence =  fdesCurrentFile.featureDescription();


    assert( specIdResult.SpectrumIdentificationItem().size() > 0 );

    addFeatureNameWithEmptyDescription( fd_sequence,"deltLCn");
    addFeatureNameWithEmptyDescription( fd_sequence,"deltCn");
    //      addFeatureNameWithEmptyDescription( fd_sequence,"Xcorr" );
    // addFeatureNameWithEmptyDescription( fd_sequence,"Sp" );
    addFeatureNameWithEmptyDescription( fd_sequence,"IonFrac" );
    addFeatureNameWithEmptyDescription( fd_sequence,"Mass" );
    addFeatureNameWithEmptyDescription( fd_sequence,"PepLen" );

    BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_cvParamType & param, specIdResult.SpectrumIdentificationItem()[0].cvParam() )  {
      if ( param.value().present() ) {
        addFeatureNameWithEmptyDescription( fd_sequence, param.name() );
      }
    }
    BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_userParamType & param, specIdResult.SpectrumIdentificationItem()[0].userParam())  {
      if ( param.value().present() ) {
        addFeatureNameWithEmptyDescription( fd_sequence, param.name() );
      }
    }
    // assert ( fdesFirstFile.featureDescription().size() > 0 || ( fdesFirstFile.featureDescription().size() == 0 && fdesCurrentFile.featureDescription().size() != fdes )) 
    if ( fdesFirstFile.featureDescription().size() == 0 ) {
      // This is the first time this function is called. We save the current Feature descriptions
      assert ( fdesCurrentFile.featureDescription().size() > 0 ); 
      fdesFirstFile = fdesCurrentFile;
    }
    // Additional files should have the same Feature descriptions as the first file
    assert ( fdesCurrentFile.featureDescription().size() == fdesFirstFile.featureDescription().size() ); 
    bool differenceFound = false;
    for ( int i = 0; i < fdesFirstFile.featureDescription().size() ;  ++i ) {
      ::percolatorInNs::featureDescription & fdesc1 = fdesFirstFile.featureDescription()[i];
      ::percolatorInNs::featureDescription & fdesc2 = fdesCurrentFile.featureDescription()[i];

      if ( fdesc1.name() != fdesc2.name() ) { differenceFound = true; break; }  
    }
    if ( differenceFound ) 
    {
      // We create the feature description list for every file. This is just a check that these list look the same.
      fprintf(stderr,"error: The file: %s translates into a feature list that is different from the a previously created feature list ( from another file ).\n", fileName); exit(EXIT_FAILURE);
    }
    for (; doc.get () != 0 && XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan>  fss_p(0);;
      scanNumberMapType::iterator iter = scanNumberMap.find( specIdResult.id() );
      int useScanNumber;
      if ( iter == scanNumberMap.end() ) {
        scanNumberMap[ specIdResult.id() ]=*scanNumber;
        useScanNumber = *scanNumber;
        ++scanNumber;
        std::auto_ptr< ::percolatorInNs::fragSpectrumScan> tmp_p( new ::percolatorInNs::fragSpectrumScan( useScanNumber, experimentalMassToCharge ));  ;
        fss_p = tmp_p;
      } else {
        useScanNumber = iter->second;
        fss_p = database.getFSS( useScanNumber );
        assert(fss_p.get());
        assert( fss_p->experimentalMassToCharge().get() == experimentalMassToCharge );
      }
      //	std::auto_ptr< ::percolatorInNs::fragSpectrumScan > fss_p
      //	::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_sequence & psm_sequence = fss_p->peptideSpectrumMatch();
      BOOST_FOREACH( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem() )  {
        createPSM(item, peptideMap, minCharge, maxCharge, experimentalMassToCharge, args_info, isDecoy, fdesFirstFile , fss_p->peptideSpectrumMatch());
      }
      database.putFSS( *fss_p );
    }
    peptideMapType::iterator iter;
    // peptideMap not needed anymore. Let us deallocate.
    for(iter = peptideMap.begin(); iter != peptideMap.end(); ++iter) { delete iter->second; iter->second=0; }
  }
  catch (const xercesc::DOMException& e)
  {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;  
    XMLString::release(&tmpStr);
    ret = 1;
  }
  catch (const xml_schema::exception& e)
  {
    cerr << e << endl;
    ret = 1;
  }
  catch (const ios_base::failure&)
  {
    cerr << "io failure" << endl;
    ret = 1;
  }
  return ret;
}

std::string greeter()
{
  ostringstream oss;
  oss << "mzidentml2pin version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2010 Lukas Käll. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukask@cbr.su.se) in the" << endl;
  oss << "Department of Biochemistry and Biophysics at the Stockholm University."
      << endl;
  return oss.str();
  
}

std::string extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << "on " << host << endl;
  return oss.str();
}

std::vector<std::string> parseFileNames(const std::string &str)
{
//   char * pch;
//   std::vector<std::string> tokens;
//   pch = strtok (str," ,|-");
//   while (pch != NULL)
//   {
//     list.insert(pch);
//     pch = strtok (NULL, " ,|-");
//   }
  
  istringstream iss(str);
  std::vector<string> tokens;
  copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(tokens));

  
  return tokens;
}

bool ParseOptions(int argc, char **argv, input_options &args_info)
{
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   mzidentml2pin [options] -t target.sqt -d decoy.sqt" << endl << endl;
  intro << "Target.sqt is the target sqt-file, and decoy.sqt is" << endl;
  intro << "the decoy sqt-file. Small data sets may be merged by replace the sqt-files with" << endl;
  intro << "meta files. Meta files are text files containing the paths of sqt-files, one" << endl;
  intro << "path per line. For successful result, the different runs should be generated" << endl;
  intro << "under similar condition." << endl;

  // init
  CommandLineParser cmd(intro.str());
  
  cmd.defineOption("t",
      "target-file",
      "A target file in mzIdentML (sequest) format",
      "filename");
  cmd.defineOption("d",
      "decoy-file",
      "A decoy file in mzIdentML (sequest) format",
      "filename");
  cmd.defineOption("e",
      "verbose",
      "Type of enzyme \"no_enzyme\",\"elastase\",\"chymotrypsin\",\"trypsin\" default=\"trypsin\"",
      "",
      "trypsin");
  cmd.defineOption("w",
      "tmp-file-for-indermediate-results",
      "tmp file for indermediate results default = /tmp/convertsequest-tmpfile.tcb",
      "filename",
      "/tmp/convertsequest-tmpfile.tcb");
  cmd.defineOption("b",
      "ptm",
      "Calculate feature for number of post-translational modifications",
      "",
      TRUE_IF_SET);
  cmd.defineOption("a",
      "aa-freq",
      "Calculate amino acid frequency features",
      "",
      TRUE_IF_SET);
  cmd.defineOption("N",
      "pngasef",
      "Calculate feature based on N-linked glycosylation pattern resulting from a PNGaseF treatment. (N[*].[ST])",
      "",
      TRUE_IF_SET);
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  
  cmd.parseArgs(argc, argv);
  
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (VERB > 0) {
    cerr << extendedGreeter();
  }
  
  args_info.tmp_file_for_indermediate_results_arg = "/tmp/convertsequest-tmpfile.tcb";
  
  if (cmd.optionSet("t")) {
    //TODO WATCH OUT IT COULD BE A LIST
    args_info.target_file_arg = parseFileNames(cmd.options["t"]);
    args_info.target_file_given = true;
  }
  if (cmd.optionSet("d")) {
    //TODO WATCH OUT IT COULD BE A LIST
    args_info.decoy_file_arg = parseFileNames(cmd.options["d"]);
    args_info.decoy_file_given = true;
  }
  if (cmd.optionSet("e")) {
    //TODO WATCH OUT IT IS A ENUM

    if( cmd.options["e"] == "no enzyme") 
      args_info.enzyme_type_arg = enzyme_type_arg_no_enzyme; 
    else if( cmd.options["e"] == "elastase") 
      args_info.enzyme_type_arg = enzyme_type_arg_elastase; 
    else if( cmd.options["e"] == "chymotrypsin")
      args_info.enzyme_type_arg = enzyme_type_arg_chymotrypsin;
    else if( cmd.options["e"] == "trypsin") 
      args_info.enzyme_type_arg = enzyme_type_arg_trypsin;
    else  
      args_info.enzyme_type_arg = enzyme_type_arg_trypsin; 

  }
  if (cmd.optionSet("w")) {
    args_info.tmp_file_for_indermediate_results_arg = cmd.options["w"];
  }
  if (cmd.optionSet("b")) {
    args_info.ptm_flag = true;
  }
  if (cmd.optionSet("a")) {
    args_info.aa_freq_flag = true;
  }
  if (cmd.optionSet("N")) {
    args_info.pngasef_flag = true;
  }
  
  return true;
}




int
main (int argc, char* argv[])
{
  /* Initialize command options parser */
  struct input_options args_info;
  if(!ParseOptions(argc, argv, args_info)){
    exit(EXIT_FAILURE);
  }
  xercesc::XMLPlatformUtils::Initialize ();

  std::vector< MinMaxStruct > vec;
  for (int i = 0; i < args_info.target_file_given; ++i) {
    getMinAndMaxCharge(args_info.target_file_arg[i].c_str(),vec);
  }
  for (int i = 0; i < args_info.decoy_file_given; ++i) {
    getMinAndMaxCharge(args_info.decoy_file_arg[i].c_str(),vec);
  }
  int minCharge;
  int maxCharge;
  BOOST_FOREACH( const MinMaxStruct & minMax, vec )  {
    minCharge = std::min(minMax.min, minCharge);
    maxCharge = std::max(minMax.max, maxCharge);
  } 
  FragSpectrumScanDatabase database;
  database.init(args_info.tmp_file_for_indermediate_results_arg);
  int ret=0;
  int scanNumber=0;
  scanNumberMapType scanNumberMap;
  std::auto_ptr<percolatorInNs::featureDescriptions> fdesFirstFile_p ( new ::percolatorInNs::featureDescriptions());
  for (int i = 0; i < args_info.target_file_given; ++i) {
    printf ("passed target file: %s\n", args_info.target_file_arg[i].c_str());
    ret = loadFromTargetOrDecoyFile(args_info.target_file_arg[i].c_str(), args_info, minCharge, maxCharge, false /* isDecoy */, *fdesFirstFile_p, database, &scanNumber,  scanNumberMap  );
    if (! ret ) { fprintf(stderr,"error: failed to read/load/parse file:%s",args_info.target_file_arg[i].c_str()); break; }
  }
  for (int i = 0; ret == 0 && i < args_info.decoy_file_given; ++i) {
    loadFromTargetOrDecoyFile(args_info.decoy_file_arg[i].c_str(), args_info, minCharge, maxCharge, true /* isDecoy */, *fdesFirstFile_p, database, &scanNumber,  scanNumberMap   );
    if (! ret ) { fprintf(stderr,"error: failed to read/load/parse file:%s",args_info.decoy_file_arg[i].c_str()); }
  }
  std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  std::cout << "<experiment  xmlns=\"" << PERCOLATOR_IN_NAMESPACE <<  "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\""  << PERCOLATOR_IN_NAMESPACE <<  " file:///scratch/e/nypercol/percolator/src/percolator-xml.xsd\">" << std::endl;

  std::string enzymeStr;
  switch(args_info.enzyme_type_arg) { 
  case enzyme_type_arg_no_enzyme: enzymeStr = "no enzyme"; break;
  case enzyme_type_arg_elastase: enzymeStr = "elastase"; break;
  case enzyme_type_arg_chymotrypsin: enzymeStr = "chymotrypsin"; break;
  case enzyme_type_arg_trypsin: enzymeStr = "trypsin"; break;
  default: break;
  }
  std::cout << "   <enzyme>" << enzymeStr << "</enzyme>" << std::endl;
  serializer ser;
  ser.start (std::cout);
  ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",  *fdesFirstFile_p );
  database.print(ser);
  std::cout << std::endl << "</experiment>" << std::endl;
  xercesc::XMLPlatformUtils::Terminate ();
  return ret;
}
