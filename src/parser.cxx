// file      : examples/cxx/tree/streaming/parser.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#include <iostream>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/dom/DOM.hpp>
#if _XERCES_VERSION >= 30000
#  include <xercesc/dom/impl/DOMTextImpl.hpp>
#endif
#include <xsd/cxx/auto-array.hxx>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/xml/sax/bits/error-handler-proxy.hxx>
#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include "parser.hxx"
#include "Globals.h"

using namespace std;
using namespace xercesc;
namespace xml = xsd::cxx::xml;
namespace tree = xsd::cxx::tree;

class parser_impl: public DefaultHandler
{
public:
  parser_impl ();

  xml::dom::auto_ptr<DOMDocument>
  start (istream& is, const string& id, bool val, string schemaDefinition,
      string schema_major, string schema_minor);

  xml::dom::auto_ptr<DOMDocument>
  next ();

  // SAX event handlers.
  //
private:
  virtual void
  startElement (const XMLCh* const uri,
                const XMLCh* const lname,
                const XMLCh* const qname,
                const Attributes& attributes);

  virtual void
  endElement (const XMLCh* const uri,
              const XMLCh* const lname,
              const XMLCh* const qname);

  virtual void
  characters (const XMLCh* const s,
#if _XERCES_VERSION >= 30000
              const XMLSize_t length
#else
              const unsigned int length
#endif
  );

private:
  // SAX parser.
  //
  bool clean_;
  auto_ptr<SAX2XMLReader> parser_;
  XMLPScanToken token_;
  tree::error_handler<char> error_handler_;
  xml::sax::bits::error_handler_proxy<char> error_proxy_;
  auto_ptr<xml::sax::std_input_source> isrc_;

  size_t depth_;
  size_t docDepth_;

  // DOM document being built.
  //
  DOMImplementation& dom_impl_;
  xml::dom::auto_ptr<DOMDocument> doc_;
  DOMElement* cur_;
};

const XMLCh ls[] = {chLatin_L, chLatin_S, chNull};

parser_impl::
parser_impl ()
    : clean_ (true),
      parser_ (XMLReaderFactory::createXMLReader ()),
      error_proxy_ (error_handler_),
      dom_impl_ (*DOMImplementationRegistry::getDOMImplementation (ls))
{
  parser_->setFeature (XMLUni::fgSAX2CoreNameSpaces, true);
  parser_->setFeature (XMLUni::fgSAX2CoreNameSpacePrefixes, true);
  parser_->setFeature (XMLUni::fgXercesValidationErrorAsFatal, true);
  parser_->setFeature (XMLUni::fgXercesSchemaFullChecking, false);

  // Xerces-C++ 3.1.0 is the first version with working multi import
  // support. It also allows us to disable buffering in the parser
  // so that the date is parsed and returned as soon as it is
  // available.
  //
#if _XERCES_VERSION >= 30100
  parser_->setFeature (XMLUni::fgXercesHandleMultipleImports, true);

  XMLSize_t lwm = 0;
  parser_->setProperty (XMLUni::fgXercesLowWaterMark, &lwm);
#endif

  parser_->setErrorHandler (&error_proxy_);
  parser_->setContentHandler (this);
}

xml::dom::auto_ptr<DOMDocument> parser_impl::
start (istream& is, const string& id, bool val, string schemaDefinition,
    string schema_major, string schema_minor)
{
  // Reset our state.
  //
  depth_ = 0;
  doc_.reset ();
  error_handler_.reset ();

  if (!clean_)
    parser_->parseReset (token_);
  else
    clean_ = false;

  isrc_.reset (new xml::sax::std_input_source (is, id));

  parser_->setFeature (XMLUni::fgSAX2CoreValidation, val);
  parser_->setFeature (XMLUni::fgXercesSchema, val);

  // if local copy of the schema is available, use it for validation...
  ifstream inp(schemaDefinition.c_str());
  if(inp && val){
    string schemaNamespace = "http://per-colator.com/percolator_in/";
    schemaNamespace.append(schema_major);
    schemaNamespace.append(schema_minor);
    schemaNamespace.append(" ");
    string schemaLocation = schemaNamespace.append(schemaDefinition);
    XMLCh* propertyValue = XMLString::transcode(schemaLocation.c_str());
    ArrayJanitor<XMLCh> janValue(propertyValue);
    parser_->
    setProperty(XMLUni::fgXercesSchemaExternalSchemaLocation, propertyValue);
  }
  // ... otherwise send a warning
  else {
    cerr << "WARNING: no valid local xml schema is available to validate the input.\n";
    cerr << "If further errors are encountered, please reinstall Percolator.\n";
  }
  
  // Start parsing. The first document that we return is a "carcase"
  // of the complete document. That is, the root element with all the
  // attributes but without any content.
  //
  bool r (parser_->parseFirst (*isrc_, token_));
  error_handler_.throw_if_failed<tree::parsing<char> > ();

  while (r && depth_ == 0)
  {
    r = parser_->parseNext (token_);
    error_handler_.throw_if_failed<tree::parsing<char> > ();
  }

  if (!r)
    return xml::dom::auto_ptr<DOMDocument> (0);

  return doc_;
}

xml::dom::auto_ptr<DOMDocument> parser_impl::
next ( )
{
  assert(doc_.get () == 0);

  //maybe remove?
  if (depth_ ==  0 )
    return xml::dom::auto_ptr<DOMDocument> (0);

  bool r (true);

  while ( r && doc_.get () == 0)
  {
    r = parser_->parseNext (token_);
    error_handler_.throw_if_failed<tree::parsing<char> > ();
  }
  if (!r)
    return xml::dom::auto_ptr<DOMDocument> (0);

  while (r && depth_ !=  docDepth_ )
  {
    r = parser_->parseNext (token_);
    error_handler_.throw_if_failed<tree::parsing<char> > ();
  }

  if (!r)
    return xml::dom::auto_ptr<DOMDocument> (0);

  return doc_;
}

void parser_impl::
startElement (const XMLCh* const uri,
              const XMLCh* const /*lname*/,
              const XMLCh* const qname,
              const Attributes& attr)
{
  if ( depth_ == 0 || depth_ == 1 || ( depth_== 4 && XMLString::equals (qname, spectrumIdentificationResultStr) ))
  {
    doc_.reset (dom_impl_.createDocument (uri, qname, 0));
    cur_ = doc_->getDocumentElement ();
    docDepth_ = depth_;
  }
  else
  {
    DOMElement* e = doc_->createElementNS (uri, qname);
    cur_->appendChild (e);
    cur_ = e;
  }

  // Set attributes.
  //
#if _XERCES_VERSION >= 30000
  for (XMLSize_t i (0), end (attr.getLength()); i < end; ++i)
#else
  for (unsigned int i (0), end (attr.getLength()); i < end; ++i)
#endif
  {
    const XMLCh* qn (attr.getQName (i));
    const XMLCh* ns (attr.getURI (i));

    // When SAX2 reports the xmlns attribute, it does not include
    // the proper attribute namespace. So we have to detect and
    // handle this case.
    //

    if (XMLString::equals (qn, XMLUni::fgXMLNSString))
       ns = XMLUni::fgXMLNSURIName;

    cur_->setAttributeNS (ns, qn, attr.getValue (i));
  }
  depth_++;
}

void parser_impl::
endElement (const XMLCh* const /*uri*/,
            const XMLCh* const /*lname*/,
            const XMLCh* const /*qname*/)
{
  // We have an element parent only on depth 2 or greater.
  //

  //  doc_.reset ();
  --depth_;

  if (doc_.get () != 0 && depth_ >= docDepth_ )
    cur_ = static_cast<DOMElement*> (cur_->getParentNode ());
}



#if _XERCES_VERSION >= 30000
void parser_impl::
characters (const XMLCh* const s, const XMLSize_t length)
{
  const XMLCh empty[] = {chNull};

  // Ignore text content (presumably whitespaces) in the root element.
  //
  if (depth_ > 1)
  {
    DOMText* t = doc_->createTextNode (empty);
    static_cast<DOMTextImpl*> (t)->appendData (s, length);
    cur_->appendChild (t);
  }
}
#else
void parser_impl::
characters (const XMLCh* const s, const unsigned int length)
{
  // Ignore text content (presumably whitespaces) in the root element.
  //
  if (depth_ > 1)
  {
    // For Xerces-C++ 2-series we have to make copy.
    //
    xsd::cxx::auto_array<XMLCh> tmp (new XMLCh[length + 1]);
    XMLString::copyNString (tmp.get (), s, length);
    cur_->appendChild (doc_->createTextNode (tmp.get ()));
  }
}
#endif


//
// parser
//

parser::
~parser ()
{
}

parser::
parser ()
    : impl_ (new parser_impl)
{
}

xml::dom::auto_ptr<DOMDocument> parser::
start (istream& is, const string& id, bool val, string schemaDefinition,
    string schema_major, string schema_minor)
{
  return impl_->start(is, id, val, schemaDefinition, schema_major, schema_minor);
}

xml::dom::auto_ptr<DOMDocument> parser::
next ()
{
  return impl_->next ();
}
