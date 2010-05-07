// file      : examples/cxx/tree/streaming/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain


// Erik Sjolund modified this file to handle the percolator xml format
// The original file comes from the streaming example in xsd-3.3.0 ( http://codesynthesis.com/download/xsd/3.3/ )
// The file examples/cxx/tree/streaming/parse.cxx could be used without changes.

#include <iostream>
#include <fstream>
#include <numeric>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <boost/foreach.hpp>

#include "parser.hxx"
#include "config.h"
#include "percolator-xml.hxx"

using namespace std;
using namespace xercesc;

void printCalibration(percolatorInNs::calibration & calibration ) {
  std::cout << "---------calibration-----------" << std::endl;
  BOOST_FOREACH( percolatorInNs::calibration_parameter par, calibration.calibration_parameter() )  {
    std::cout << "name=" << par.name() << " value=" << par.value() << std::endl;
  }  
  std::cout << "mass_type=" << calibration.mass_type() << std::endl;
}

void exampleUsage(percolatorInNs::calibration & calibration ) {
  if ( percolatorInNs::mass_type::monoisotopic == calibration.mass_type() ) {
    std::cout << "percolatorInNs::mass_type::monoisotopic == calibration.mass_type() is true" << std::endl;  
  }
  if ( percolatorInNs::mass_type::average == calibration.mass_type() ) {
    std::cout << "percolatorInNs::mass_type::average == calibration.mass_type() is true" << std::endl;  
  }
}

void exampleUsage(percolatorInNs::peptide_spectrum_match & psm ) {
  percolatorInNs::features::feature_sequence & v =  psm.features().feature();
  cout << "Sum of the features=" <<  std::accumulate( v.begin(), v.end(), 0.0 )  << std::endl; 
}

void printFeatureDescriptions(percolatorInNs::feature_descriptions & feature_descriptions ) {
  std::cout << "--------- feature_descriptions -----------" << std::endl;
  BOOST_FOREACH( percolatorInNs::feature_descriptions::feature_description_type fdes, feature_descriptions.feature_description()  )  {
    std::cout << fdes << std::endl;
  }
}

void printFragSpectrumScan(percolatorInNs::frag_spectrum_scan &fss) {
  std::cout << "--------- frag_spectrum_scan -----------" << std::endl;
  std::cout << " num=" << fss.num() << " observed mass charge=" << fss.observed().mass_charge() << std::endl;
  if (fss.ion_current().present()) { 
    std::cout << " ion current=" << fss.ion_current().get().val() << std::endl;
  }
  BOOST_FOREACH( percolatorInNs::peptide_spectrum_match psm, fss.peptide_spectrum_match() )  {
    exampleUsage(psm);
    std::cout << " charge=" << psm.charge() << " type=" << psm.type() << " id=" << psm.id() << std::endl;
    BOOST_FOREACH( percolatorInNs::features::feature_type feature, psm.features().feature() )  {
           std::cout << "  feature=" << feature << std::endl;
    }
    std::cout << " calculated mass charge=" << psm.calculated().mass_charge() << " peptide=" << psm.peptide() << std::endl;
  }
  // retention_time is optional so we first need to check its presence
  if (fss.observed().retention_time().present ()) { 
    std::cout << "retention_time=" <<  fss.observed().retention_time().get() << std::endl;
  }
} 

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

  try
  {
    namespace xml = xsd::cxx::xml;
    ifstream ifs;
    ifs.exceptions (ifstream::badbit | ifstream::failbit);
    ifs.open (argv[1]);

    parser p;

    xml_schema::dom::auto_ptr<DOMDocument> doc (p.start (ifs, argv[1], true));


    // The num_features attribute is not in one of the sub trees so it need some special treatment
    DOMAttr* num_features_attr ( doc->getDocumentElement ()->getAttributeNode ( xml::string ("num_features").c_str ()));
    percolatorInNs::experiment::num_features_type num_f ( percolatorInNs::experiment::num_features_traits::create (*num_features_attr, 0, 0));
    std::cout << "num_features=" << num_f << std::endl;

    doc = p.next ();

    // The enzyme element is a subelement but CodeSynthesis Xsd does not generate a class for it. (I am trying to find a command line option that overrides this decision)
    // As for now special treatment is needed:
    char * value = XMLString::transcode(   doc->getDocumentElement()->getTextContent()  );
    std::cout << "enzyme=" <<  value << std::endl;   
    XMLString::release(&value);

    doc = p.next ();

    static const XMLCh calibrationStr[] = { chLatin_c, chLatin_a, chLatin_l, chLatin_i, chLatin_b,chLatin_r, chLatin_a, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };
    if (XMLString::equals( calibrationStr, doc->getDocumentElement ()->getTagName())) {  
      percolatorInNs::calibration calibration(*doc->getDocumentElement ());
      printCalibration(calibration);
      exampleUsage(calibration);
      doc = p.next ();
    };

    percolatorInNs::feature_descriptions feature_descriptions(*doc->getDocumentElement ());
    printFeatureDescriptions(feature_descriptions);

    for (doc = p.next (); doc.get () != 0; doc = p.next ())
    {
      percolatorInNs::frag_spectrum_scan frag_spectrum_scan(*doc->getDocumentElement ());
      printFragSpectrumScan(frag_spectrum_scan);
    }
  }
  catch (const xercesc_3_1::DOMException& e)
  {
    char * tmpStr = XMLString::transcode(e.getMessage());
    std::cerr << "catch  xercesc_3_1::DOMException=" << tmpStr << std::endl;  
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
