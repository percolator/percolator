#include "MzidentmlReader.h"
#include "DataSet.h"

static const XMLCh sequenceCollectionStr[] = { chLatin_S, chLatin_e, chLatin_q, chLatin_u, chLatin_e,chLatin_n, 
						  chLatin_c, chLatin_e, chLatin_C, chLatin_o, chLatin_l, chLatin_l, chLatin_e, 
						  chLatin_c, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };
						  
						  

MzidentmlReader::MzidentmlReader(ParseOptions po):Reader(po)
{

}

MzidentmlReader::~MzidentmlReader()
{

}

bool MzidentmlReader::checkValidity(string file)
{
  bool ismeta = false;
  std::ifstream fileIn(file.c_str(), std::ios::in);
  if (!fileIn) {
    std::cerr << "Could not open file " << file << std::endl;
    exit(-1);
  }
  std::string line;
  if (!getline(fileIn, line)) {
    std::cerr << "Could not read file " << file << std::endl;
    exit(-1);
  }
  fileIn.close();
  if (line.size() > 1 && line[0]=='<' && line[1]=='?') {
    //NOTE I should check deeper and validate that the file is mzidentml
    if (line.find("xml") == std::string::npos) {
      std::cerr << "file is not xml format " << file << std::endl;
      exit(-1);
    }
  }
  else
  {
    ismeta = true;
  }
  return ismeta;
}


void MzidentmlReader::addFeatureDescriptions(bool doEnzyme, const string& aaAlphabet, std::string fn)
{

  ifstream ifs;
  ifs.exceptions (ifstream::badbit | ifstream::failbit);
  try
  {
    ifs.open (fn.c_str());
    parser p;
    string schemaDefinition = PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, fn, true, schemaDefinition,
      schema_major, schema_minor));
    for (doc = p.next (); doc.get () != 0 && 
      !XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() );
      doc = p.next ()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }
    ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
    
    push_backFeatureDescription("deltLCn");
    push_backFeatureDescription("deltCn");
    //push_backFeatureDescription("Xcorr");
    //push_backFeatureDescription("Sp");
    push_backFeatureDescription("IonFrac");
    push_backFeatureDescription("Mass" );
    push_backFeatureDescription("PepLen");


    BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_cvParamType & param, specIdResult.SpectrumIdentificationItem()[0].cvParam() )  {
      if ( param.value().present() ) {
        push_backFeatureDescription(param.name().c_str());
      }
    }
    BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_userParamType & param, specIdResult.SpectrumIdentificationItem()[0].userParam())  {
      if ( param.value().present() ) {
        push_backFeatureDescription(param.name().c_str());
      }
    }
    
    
    //NOTE are these added already??
    /*if (doEnzyme) {
      push_backFeatureDescription("enzN");
      push_backFeatureDescription("enzC");
      push_backFeatureDescription("enzInt");
    }

    if (po.calcPTMs) {
      push_backFeatureDescription("ptm");
    }
    if (po.pngasef) {
      push_backFeatureDescription("PNGaseF");
    }
    if (!aaAlphabet.empty()) {
      for (std::string::const_iterator it = aaAlphabet.begin(); it
      != aaAlphabet.end(); it++)
	push_backFeatureDescription(*it + "-Freq");
    }*/

  }catch (ifstream::failure e) {
    cerr << "Exception opening/reading file :" << fn <<endl;
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
  
  return; 
}


void MzidentmlReader::getMaxMinCharge(string fn)
{

  bool foundFirstChargeState = false;
  ifstream ifs;
  ifs.exceptions (ifstream::badbit | ifstream::failbit);
  try
  {
    ifs.open (fn.c_str());
    parser p;
    string schemaDefinition = PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, fn, true, schemaDefinition,
      schema_major, schema_minor));
    for (doc = p.next (); doc.get () != 0 && 
      !XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() );
      doc = p.next ()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }
    ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());

    int scanNumber = 0;
    for (; doc.get () != 0 && XMLString::equals( spectrumIdentificationResultStr, 
           doc->getDocumentElement ()->getTagName() ); doc = p.next ()) 
    {
      ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
      
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();
      
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan > fss_p( new ::percolatorInNs::fragSpectrumScan( scanNumber, experimentalMassToCharge )); 
      
      scanNumber++;
      BOOST_FOREACH( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem() )  {
	if ( ! foundFirstChargeState ) {
	  minCharge = item.chargeState();
	  minCharge = item.chargeState();
	  foundFirstChargeState = true;
	}
	minCharge = std::min(item.chargeState(),minCharge);
	minCharge = std::max(item.chargeState(),minCharge);
      }
    }
  }catch (ifstream::failure e) {
    cerr << "Exception opening/reading file :" << fn <<endl;
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
  return;
}

void MzidentmlReader::read(const std::string fn, bool isDecoy, boost::shared_ptr<FragSpectrumScanDatabase> database) 
{
  namespace xml = xsd::cxx::xml;
  int scanNumber=0;
  scanNumberMapType scanNumberMap;
  
  try
  {
    ifstream ifs;
    ifs.exceptions (ifstream::badbit | ifstream::failbit);
    ifs.open (fn.c_str());
    parser p;
    string schemaDefinition = PIN_SCHEMA_LOCATION + string("percolator_in.xsd");
    string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
    string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, fn.c_str(), true, schemaDefinition,
        schema_major, schema_minor));
    
    while (doc.get () != 0 && ! XMLString::equals( sequenceCollectionStr, doc->getDocumentElement ()->getTagName())) {
      doc = p.next ();
      // Let's skip some sub trees that we are not interested, e.g. AuditCollection
    }
    assert(doc.get());
    mzIdentML_ns::SequenceCollectionType sequenceCollection(*doc->getDocumentElement ());

    peptideMapType peptideMap;

    BOOST_FOREACH( const mzIdentML_ns::SequenceCollectionType::Peptide_type &peptide, sequenceCollection.Peptide() )  {
      assert( peptideMap.find( peptide.id() ) == peptideMap.end() ); // The peptide refs should be unique.
      mzIdentML_ns::SequenceCollectionType::Peptide_type *pept = new mzIdentML_ns::SequenceCollectionType::Peptide_type( peptide);
      assert(pept);
      peptideMap.insert( std::make_pair(peptide.id(), pept ) ) ;
    }
    
    for (doc = p.next (); doc.get () != 0 && !XMLString::equals( spectrumIdentificationResultStr, 
      doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      // Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
    }
    
    ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
    assert( specIdResult.SpectrumIdentificationItem().size() > 0 );

    for (; doc.get () != 0 && XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
      ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
      assert(specIdResult.SpectrumIdentificationItem().size() > 0);
      ::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan>  fss_p(0);;
      scanNumberMapType::iterator iter = scanNumberMap.find( specIdResult.id() );
      int useScanNumber;
      if ( iter == scanNumberMap.end() ) {
        scanNumberMap[ specIdResult.id() ]=scanNumber;
        useScanNumber = scanNumber;
        ++scanNumber;
        std::auto_ptr< ::percolatorInNs::fragSpectrumScan> tmp_p( new ::percolatorInNs::fragSpectrumScan( useScanNumber, experimentalMassToCharge ));  ;
        fss_p = tmp_p;
      } else {
        useScanNumber = iter->second;
        fss_p = database->getFSS( useScanNumber );
        assert(fss_p.get());
        assert( fss_p->experimentalMassToCharge().get() == experimentalMassToCharge );
      }

      BOOST_FOREACH( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, specIdResult.SpectrumIdentificationItem() )  {
        createPSM(item, peptideMap, experimentalMassToCharge, isDecoy, fss_p->peptideSpectrumMatch());
      }
      database->putFSS( *fss_p );
    }
    
    peptideMapType::iterator iter;
    for(iter = peptideMap.begin(); iter != peptideMap.end(); ++iter) 
    { 
      delete iter->second; iter->second=0;
    }
  }
  catch (const xercesc::DOMException& e)
  {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;  
    XMLString::release(&tmpStr);
    exit(-1);
  }
  catch (const xml_schema::exception& e)
  {
    cerr << e << endl;
    exit(-1);
  }
  catch (const ios_base::failure&)
  {
    cerr << "io failure" << endl;
    exit(-1);
  }
}

void MzidentmlReader::createPSM(const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, 
		peptideMapType & peptideMap,::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type 
		experimentalMassToCharge,bool isDecoy,::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_sequence & psm_sequence) 
{
  assert( item.experimentalMassToCharge() == experimentalMassToCharge );
  assert( item.Peptide_ref().present() );
  assert( peptideMap.find( item.Peptide_ref().get() ) != peptideMap.end() );
  std::string peptideSeq =  peptideMap[item.Peptide_ref().get()]->peptideSequence();
  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideSeq ) );
  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
  double rank = item.rank();

  if ( ! item.calculatedMassToCharge().present() ) 
  { 
    std::cerr << "error: calculatedMassToCharge attribute is needed for percolator" << std::endl; 
    exit(EXIT_FAILURE); 
  }
  
  f_seq.push_back( 0.0 ); // delt5Cn (leave until last M line)
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

  if ( Enzyme::getEnzymeType() != Enzyme::NO_ENZYME ) {
    f_seq.push_back( Enzyme::isEnzymatic(peptideSeqWithFlanks.at(0),peptideSeqWithFlanks.at(2)) ? 1.0
		     : 0.0);
    f_seq.push_back( 
		    Enzyme::isEnzymatic(peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 3),
					peptideSeqWithFlanks.at(peptideSeqWithFlanks.size() - 1))
		    ? 1.0
		    : 0.0);
    f_seq.push_back( (double)Enzyme::countEnzymatic(peptideSeq) );
  }
  f_seq.push_back( dM ); // obs - calc mass
  f_seq.push_back( (dM < 0 ? -dM : dM)); // abs only defined for integers on some systems
  if (po.calcPTMs ) { f_seq.push_back(  DataSet::cntPTMs(peptideSeqWithFlanks)); }
  if (po.pngasef ) { f_seq.push_back( DataSet::isPngasef(peptideSeqWithFlanks, isDecoy)); }
  if (po.calcAAFrequencies ) { computeAAFrequencies(peptideSeqWithFlanks, f_seq); }
   
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
  assert( this->f_seq.featureDescription().size() == f_seq.size() );
  //is optional for mzIdentML1.0.0 but is compulsory for percolator_in ..... Is this ok?
  assert( item.calculatedMassToCharge().present() );
  ::percolatorInNs::peptideSpectrumMatch* tmp_psm =
      new ::percolatorInNs::peptideSpectrumMatch( features_p, peptide_p, item.id(), isDecoy, item.experimentalMassToCharge(), item.calculatedMassToCharge().get(), item.chargeState());
  std::auto_ptr< ::percolatorInNs::peptideSpectrumMatch > psm_p(tmp_psm);
  //NOTE I should do this in the database object
  psm_sequence.push_back(psm_p);
  return;
}