// Erik Sjolund modified this file to handle the percolator xml format
// The original file comes from the streaming example in xsd-3.3.0 ( http://codesynthesis.com/download/xsd/3.3/ )

// original file      : examples/cxx/tree/streaming/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>

#include <iostream>
#include <fstream>
#include <numeric>
#include <map>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <boost/foreach.hpp>

#include "parser.hxx"
#include "config.h"
#include "percolator_in.hxx"
#include "mzIdentML1.0.0.hxx"
#include "serializer.hxx"
          
using namespace std;
using namespace xercesc;

/* A sketchy overview of the conversion. ( Xpath is used in the explanation )

First a hash is created with 

/mzIdentML/SequenceCollection/Peptide/@id 

as key and the subtree 

/mzIdentML/SequenceCollection/Peptide 

as value.

Then each 

/mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult  

will be read into memory and translated into a 

/experiment/fragSpectrumScan

Although a 

/mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/@id 

exists, it's a string and not a number so we can't use it. Instead we increment

/experiment/fragSpectrumScan/@scan_number  

with +1 for each new /experiment/fragSpectrumScan/  ( starting with 0 )

Each
/mzIdentML/DataCollection/AnalysisData/SpectrumIdentificationList/SpectrumIdentificationResult/SpectrumIdentificationItem
translates into a 
/experiment/fragSpectrumScan/peptideSpectrumMatch


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

---------------------

The next() method of parser p, needs an explanation:

It returns:
first time just the root element
then the next subtree child of the root element or the next SpectrumIdentificationResult sub tree.
*/

int
main (int argc, char* argv[])
{
  if (argc != 2)
    {
      cerr << "usage: " << argv[0] << " percolator-in.xml" << endl;
      return 1;
    }
  int r (0);
  xercesc::XMLPlatformUtils::Initialize ();
  bool isDecoy = true;
  std::string enzyme="Mitt enzym";
 try
    {
      namespace xml = xsd::cxx::xml;
      ifstream ifs;
      ifs.exceptions (ifstream::badbit | ifstream::failbit);
      ifs.open (argv[1]);
      parser p;
      xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, argv[1], true));
      std::auto_ptr<percolatorInNs::featureDescriptions> fdes_p ( new ::percolatorInNs::featureDescriptions());
      std::auto_ptr< ::percolatorInNs::experiment > ex_p ( new ::percolatorInNs::experiment( "mitt enzym" , fdes_p ));
      static const XMLCh sequenceCollectionStr[] = { chLatin_S, chLatin_e, chLatin_q, chLatin_u, chLatin_e,chLatin_n, chLatin_c, chLatin_e, chLatin_C, chLatin_o, chLatin_l, chLatin_l, chLatin_e, chLatin_c, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };
      std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
      std::cout << "<experiment  xmlns=\"" << PERCOLATOR_IN_NAMESPACE <<  "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\""  << PERCOLATOR_IN_NAMESPACE <<  " file:///scratch/e/nypercol/percolator/src/percolator-xml.xsd\">" << std::endl;
      std::cout << "   <enzyme>" << enzyme << "</enzyme>" << std::endl;
      while (doc.get () != 0 && ! XMLString::equals( sequenceCollectionStr, doc->getDocumentElement ()->getTagName())) {
	doc = p.next ();
	// Let's skip some sub trees that we are not interested, e.g. AuditCollection
      }
      assert(doc.get());
      mzIdentML_ns::SequenceCollectionType sequenceCollection(*doc->getDocumentElement ());
      // instead of std::map it should actually be a unordered_map. Let us wait until c++0x finalize.
      typedef map<std::string, mzIdentML_ns::SequenceCollectionType::Peptide_type *> peptideMapType;
      peptideMapType peptideMap;
      BOOST_FOREACH( const mzIdentML_ns::SequenceCollectionType::Peptide_type peptide, sequenceCollection.Peptide() )  {
        assert( peptideMap.find( peptide.id() ) == peptideMap.end() ); // The peptide refs should be unique.
	mzIdentML_ns::SequenceCollectionType::Peptide_type *pept = new mzIdentML_ns::SequenceCollectionType::Peptide_type( peptide);
        assert(pept);
        peptideMap.insert( std::make_pair(peptide.id(), pept ) ) ;
      }
      for (doc = p.next (); doc.get () != 0 && !XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
	// Let's skip some sub trees that we are not interested, e.g. AnalysisCollection
      }
      serializer ser;
      ser.start (std::cout);
      ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
      percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence =  ex_p->featureDescriptions().featureDescription();
      assert( specIdResult.SpectrumIdentificationItem().size() > 0 );
      BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_cvParamType cv, specIdResult.SpectrumIdentificationItem()[0].cvParam() )  {
        if ( cv.value().present() ) {
	  std::auto_ptr< ::percolatorInNs::featureDescription > f_p( new ::percolatorInNs::featureDescription( cv.name() )); 
	  fd_sequence.push_back( f_p );
	}
      }
      ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions", ex_p->featureDescriptions() );
      int scanNumber = 0;
      for (; doc.get () != 0 && XMLString::equals( spectrumIdentificationResultStr, doc->getDocumentElement ()->getTagName() ); doc = p.next ()) {
	::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationResultType specIdResult(*doc->getDocumentElement ());
	assert(specIdResult.SpectrumIdentificationItem().size() > 0);
	::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type experimentalMassToCharge = specIdResult.SpectrumIdentificationItem()[0].experimentalMassToCharge();
	std::auto_ptr< ::percolatorInNs::fragSpectrumScan > fss_p( new ::percolatorInNs::fragSpectrumScan( scanNumber, experimentalMassToCharge )); 
	scanNumber++;
	BOOST_FOREACH( const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType item, specIdResult.SpectrumIdentificationItem() )  {
          // It is strange but mzIdentML has "experimentalMassToCharge" on the PSM-level in the XML tree. This leads to a lot of redundant information.
          // Let us check that our assumption ( experimentalMassToCharge is constant in the subtree under SpectrumIdentificationResult ) is really valid with an assert()
          assert( item.experimentalMassToCharge() == experimentalMassToCharge );
       	  assert( item.Peptide_ref().present() );
	  assert( peptideMap.find( item.Peptide_ref().get() ) != peptideMap.end() );
	  std::auto_ptr< percolatorInNs::peptideType >  peptide_p( new percolatorInNs::peptideType( peptideMap[item.Peptide_ref().get()]->peptideSequence() ) );
	  std::auto_ptr< percolatorInNs::features >  features_p( new percolatorInNs::features ());
	  percolatorInNs::features::feature_sequence & f_seq =  features_p->feature();
	  BOOST_FOREACH( const ::mzIdentML_ns::FuGE_Common_Ontology_cvParamType cv, item.cvParam() )  {
	    if ( cv.value().present() ) {
	      // SpectrumIdentificationItem/cvParam/@value has the datatype string, even though the values seem to be float or double. 
              // percolator_in.xsd uses the datatype double for the features/feature, so we need to convert the string.
              // Using feature_traits for the conversion from std::string to double seems to be the right way to go. Another option would have been to use "strtod()".
	      percolatorInNs::features::feature_type fe ( percolatorInNs::features::feature_traits::create(cv.value().get().c_str(),0,0,0) );
	      f_seq.push_back( fe);
	    }
	  }
	  assert( fd_sequence.size() == f_seq.size() );
	  //is optional for mzIdentML1.0.0 but is compulsory for percolator_in ..... Is this ok?
	  assert( item.calculatedMassToCharge().present() );
	  std::auto_ptr< ::percolatorInNs::peptideSpectrumMatch > psm_p( new ::percolatorInNs::peptideSpectrumMatch( features_p, peptide_p, item.id(), isDecoy, item.calculatedMassToCharge().get(), item.chargeState() )); 
	  fss_p->peptideSpectrumMatch().push_back(psm_p);
	}
        ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss_p );
      }
      std::cout << std::endl << "</experiment>" << std::endl;
      peptideMapType::iterator iter;
      // peptideMap not needed anymore. Let us deallocate.
      for(iter = peptideMap.begin(); iter != peptideMap.end(); ++iter) { delete iter->second; iter->second=0; }
   }
  catch (const xercesc_3_1::DOMException& e)
    {
      char * tmpStr = XMLString::transcode(e.getMessage());
      std::cerr << "catch xercesc_3_1::DOMException=" << tmpStr << std::endl;  
      XMLString::release(&tmpStr);
      r = 1;
    }
  catch (const xml_schema::exception& e)
    {
      cerr << e << endl;
      r = 1;
    }
  catch (const ios_base::failure&)
    {
      cerr << "io failure" << endl;
      r = 1;
    }
  xercesc::XMLPlatformUtils::Terminate ();
  return r;
}
