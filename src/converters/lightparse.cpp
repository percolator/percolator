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
// #include <rpc/types.h>
// #include <rpc/xdr.h>
// #if defined __WIN32
//  #include <xdr_api_mingw.h>
// #endif
#include "percolator_in.hxx"

using namespace std;
using namespace xercesc;

void printCalibration(const percolatorInNs::calibration & calibration ) {
  std::cout << "---------calibration-----------" << std::endl;
  BOOST_FOREACH( const percolatorInNs::calibrationParameter & par, calibration.calibrationParameter() )  {
    std::cout << "name=" << par.name() << " value=" << par.value() << std::endl;
  }  
  std::cout << "massType=" << calibration.massType() << std::endl;
}

void exampleUsage(const percolatorInNs::calibration & calibration ) {
  if ( percolatorInNs::massType::monoisotopic == calibration.massType() ) {
    std::cout << "percolatorInNs::mass_type::monoisotopic == calibration.massType() is true" << std::endl;  
  }
  if ( percolatorInNs::massType::average == calibration.massType() ) {
    std::cout << "percolatorInNs::mass_type::average == calibration.massType() is true" << std::endl;  
  }
}

void exampleUsage(const percolatorInNs::peptideSpectrumMatch & psm ) {
  const percolatorInNs::features::feature_sequence & v =  psm.features().feature();
  cout << "Sum of the features=" <<  std::accumulate( v.begin(), v.end(), 0.0 )  << std::endl; 
}

void printFeatureDescriptions(const percolatorInNs::featureDescriptions & feature_descriptions ) {
  std::cout << "--------- featureDescriptions -----------" << std::endl;
  BOOST_FOREACH( const percolatorInNs::featureDescriptions::featureDescription_type & fdes, feature_descriptions.featureDescription()  )  {
     std::cout << "  feature name=" << fdes.name()  << std::endl;
  }
}

void printFragSpectrumScan(percolatorInNs::fragSpectrumScan &fss) {
  std::cout << "--------- fragSpectrumScan -----------" << std::endl;
  std::cout << " num=" << fss.scanNumber() << " experimental mass to charge=" << fss.experimentalMassToCharge() << std::endl;
  if (fss.totalIonCurrent().present()) { 
    std::cout << " ion current=" << fss.totalIonCurrent().get() << std::endl;
  }
  BOOST_FOREACH( const percolatorInNs::peptideSpectrumMatch & psm, fss.peptideSpectrumMatch() )  {
    exampleUsage(psm);
    if ( psm.isDecoy() ) {
    std::cout << " PSM is a decoy!" << std::endl;
    }
    std::cout << " charge_state=" << psm.chargeState() << " id=" << psm.id() << std::endl;
    BOOST_FOREACH( const percolatorInNs::features::feature_type & feature, psm.features().feature() )  {
      std::cout << "  feature=" << feature << std::endl;
    }
    std::cout << " calculated mass to charge=" << psm.calculatedMassToCharge() << " peptide=" << psm.peptide().peptideSequence() << std::endl;
  }
  // retention_time is optional so we first need to check its presence
  if (fss.observedTime().present ()) { 
    std::cout << "observed time=" <<  fss.observedTime().get() << std::endl;
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



    doc = p.next();

    // The enzyme element is a subelement but CodeSynthesis Xsd does not generate a class for it. (I am trying to find a command line option that overrides this decision)
    // As for now special treatment is needed:
    char * value = XMLString::transcode(   doc->getDocumentElement()->getTextContent()  );
    std::cout << "enzyme=" <<  value << std::endl;   
    XMLString::release(&value);

    doc = p.next();

    static const XMLCh calibrationStr[] = { chLatin_c, chLatin_a, chLatin_l, chLatin_i, chLatin_b,chLatin_r, chLatin_a, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };
    if (XMLString::equals( calibrationStr, doc->getDocumentElement ()->getTagName())) {  
      percolatorInNs::calibration calibration(*doc->getDocumentElement ());
      printCalibration(calibration);
      exampleUsage(calibration);
      doc = p.next();
    };

    percolatorInNs::featureDescriptions feature_descriptions(*doc->getDocumentElement ());
    printFeatureDescriptions(feature_descriptions);

    for (doc = p.next(); doc.get () != 0; doc = p.next())
    {
      percolatorInNs::fragSpectrumScan frag_spectrum_scan(*doc->getDocumentElement ());
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
