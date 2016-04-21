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
#include "MyException.h"

#include <string>
#include <memory>   // std::auto_ptr
#include <cstddef>  // std::size_t

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XercesVersion.hpp>

#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>

#include <xercesc/validators/common/Grammar.hpp>
#include <xercesc/framework/XMLGrammarPoolImpl.hpp>

using namespace std;
using namespace xercesc;
namespace xml = xsd::cxx::xml;
namespace tree = xsd::cxx::tree;

class error_handler: public ErrorHandler {
 public:
  error_handler() : failed_ (false) {}

  bool failed() const { return failed_; }

  enum severity {s_warning, s_error, s_fatal};

  virtual void warning(const SAXParseException&);
  virtual void error(const SAXParseException&);
  virtual void fatalError(const SAXParseException&);
  virtual void resetErrors() { failed_ = false; }

  void handle(const SAXParseException&, severity);

 private:
  bool failed_;
};

class parser_impl : public DefaultHandler {
 public:
  parser_impl();

  xml::dom::auto_ptr<DOMDocument> start (istream& is, const string& id, 
    bool val, string schemaDefinition, string schema_major, string schema_minor, 
    string schemaNamespace = "http://per-colator.com/percolator_in/", 
    bool noNameSpace = false);
  xml::dom::auto_ptr<DOMDocument> next();

  // SAX event handlers.
 private:
  virtual void startElement(const XMLCh* const uri,
                            const XMLCh* const lname,
                            const XMLCh* const qname,
                            const Attributes& attributes);

  virtual void endElement(const XMLCh* const uri,
                          const XMLCh* const lname,
                          const XMLCh* const qname);

  virtual void characters(const XMLCh* const s,
#if _XERCES_VERSION >= 30000
                          const XMLSize_t length
#else
                          const unsigned int length
#endif
  );

 private:
  // SAX parser.
  bool clean_;
  XMLPScanToken token_;
  error_handler error_handler_;
  std::auto_ptr<xml::sax::std_input_source> isrc_;

  size_t depth_;
  size_t docDepth_;

  // DOM document being built.
  DOMImplementation& dom_impl_;
  xml::dom::auto_ptr<DOMDocument> doc_;
  DOMElement* cur_;
  MemoryManager* mm_;// (XMLPlatformUtils::fgMemoryManager);
  std::auto_ptr<XMLGrammarPool> gp_; // (new XMLGrammarPoolImpl (mm));
  std::auto_ptr<SAX2XMLReader> parser_; // (create_parser (gp.get ()));
};

const XMLCh ls[] = {chLatin_L, chLatin_S, chNull};

std::auto_ptr<SAX2XMLReader> create_parser(XMLGrammarPool* pool) {
  std::auto_ptr<SAX2XMLReader> parser(
    pool
    ? XMLReaderFactory::createXMLReader(XMLPlatformUtils::fgMemoryManager, pool)
    : XMLReaderFactory::createXMLReader());

  // Commonly useful configuration.
  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);
  parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, true);
  parser->setFeature(XMLUni::fgSAX2CoreValidation, true);

  // Enable validation.
  parser->setFeature(XMLUni::fgXercesSchema, true);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);
  parser->setFeature(XMLUni::fgXercesValidationErrorAsFatal, true);

  // Use the loaded grammar during parsing.
  parser->setFeature(XMLUni::fgXercesUseCachedGrammarInParse, true);

  // Don't load schemas from any other source (e.g., from XML document's
  // xsi:schemaLocation attributes).
  parser->setFeature(XMLUni::fgXercesLoadSchema, false);

  // Xerces-C++ 3.1.0 is the first version with working multi import
  // support.
#if _XERCES_VERSION >= 30100
  parser->setFeature(XMLUni::fgXercesHandleMultipleImports, true);
#endif

  return parser;
}

parser_impl::parser_impl() : clean_(true), 
    mm_(XMLPlatformUtils::fgMemoryManager), gp_(new XMLGrammarPoolImpl(mm_)),
    dom_impl_(*DOMImplementationRegistry::getDOMImplementation(ls)) {}


xml::dom::auto_ptr<DOMDocument> parser_impl::start(istream& is, 
    const string& id, bool val, string schemaDefinition, string schema_major, 
    string schema_minor,string schemaNamespace, bool noNameSpace) {
  // Reset our state.
  depth_ = 0;
  doc_.reset();
  error_handler_.resetErrors();
  if (!clean_)
    parser_->parseReset(token_);
  else
    clean_ = false;
  
  if (id.size() > 0) {
    isrc_.reset (new xml::sax::std_input_source(is, id));
  } else {
    isrc_.reset (new xml::sax::std_input_source(is));
  }

  int r(0);

  while (true) {
    int i(1);
    
    // Load the schemas into the grammar pool.
    {
      std::auto_ptr<SAX2XMLReader> parser(create_parser(gp_.get()));
      parser->setErrorHandler(&error_handler_);
      
      string s(schemaDefinition);
      size_t n(s.size());

      if (!parser->loadGrammar(s.c_str(), Grammar::SchemaGrammarType, true)) {
        cerr << s << ": error: unable to load" << endl;
        r = 1;
        break;
      }/*
      if (eh.failed()) {
        r = 1;
        break;
      }*/

      if (r != 0) break;
    }
    // Lock the grammar pool. This is necessary if we plan to use the
    // same grammar pool in multiple threads (this way we can reuse the
    // same grammar in multiple parsers). Locking the pool disallows any
    // modifications to the pool, such as an attempt by one of the threads
    // to cache additional schemas.
    gp_->lockPool();

    // Parse the XML documents.
    parser_ = create_parser(gp_.get());  
    parser_->setErrorHandler(&error_handler_);  
    parser_->setContentHandler(this);
    
    // Start parsing. The first document that we return is a "carcase"  
    // of the complete document. That is, the root element with all the
    // attributes but without any content.
    bool r (parser_->parseFirst(*isrc_, token_));

    while (r && depth_ == 0) {
      r = parser_->parseNext (token_);
    }

    if (!r)
      return xml::dom::auto_ptr<DOMDocument>(0);

    /*if (eh.failed ()) {
      r = 1;
      break;
    }*/
    break;
  }
  return doc_;
}

void error_handler::warning(const SAXParseException& e) {
  handle(e, s_warning);
}

void error_handler::error(const SAXParseException& e) {
  failed_ = true;
  handle(e, s_error);
}

void error_handler::fatalError(const SAXParseException& e) {
  failed_ = true;
  handle(e, s_fatal);
}

void error_handler::handle(const SAXParseException& e, severity s) {
  const XMLCh* xid(e.getPublicId());

  if (xid == 0)
    xid = e.getSystemId();
  
  if (xid == 0) // MT: this should not happen, but it does...
    xid = XMLString::transcode("input.xml");
  
  char* id(XMLString::transcode(xid));
  char* msg(XMLString::transcode(e.getMessage()));
  
  ostringstream temp;
  temp << "XML parser " << (s == s_warning ? "warning" : "error") 
     << " at " << id << ":" << e.getLineNumber() << ":" << e.getColumnNumber()
     << "\n  " << (s == s_warning ? "warning: " : "error: ") << msg << endl;
  XMLString::release (&id);
  XMLString::release (&msg);
  
  if (s == s_fatal || s == s_error) {
    throw MyException(temp.str());
  } else {
    std::cerr << temp.str();
  }
}


xml::dom::auto_ptr<DOMDocument> parser_impl::next() {
  assert(doc_.get() == 0);

  //maybe remove?
  if (depth_ == 0)
    return xml::dom::auto_ptr<DOMDocument>(0);

  bool r (true);

  while (r && doc_.get() == 0) {
    r = parser_->parseNext (token_);
    //error_handler_.throw_if_failed<tree::parsing<char> >();
  }
  if (!r)
    return xml::dom::auto_ptr<DOMDocument>(0);

  while (r && depth_ != docDepth_) {
    r = parser_->parseNext (token_);
    //error_handler_.throw_if_failed<tree::parsing<char> >();
  }

  if (!r)
    return xml::dom::auto_ptr<DOMDocument>(0);

  return doc_;
}

void parser_impl::startElement(const XMLCh* const uri,
    const XMLCh* const /*lname*/, const XMLCh* const qname,
    const Attributes& attr) {
  //std::cerr << depth_ << " " << XMLString::transcode(qname) << std::endl;
  if (depth_ == 0 || depth_ == 1 || 
      (depth_== 4 && XMLString::equals(qname, spectrumIdentificationResultStr))) {
    doc_.reset(dom_impl_.createDocument(uri, qname, 0));
    cur_ = doc_->getDocumentElement();
    docDepth_ = depth_;
  } else {
    DOMElement* e = doc_->createElementNS(uri, qname);
    cur_->appendChild(e);
    cur_ = e;
  }

  // Set attributes.
#if _XERCES_VERSION >= 30000
  for (XMLSize_t i(0), end(attr.getLength()); i < end; ++i)
#else
  for (unsigned int i(0), end(attr.getLength()); i < end; ++i)
#endif
  {
    const XMLCh* qn(attr.getQName(i));
    const XMLCh* ns(attr.getURI(i));

    // When SAX2 reports the xmlns attribute, it does not include
    // the proper attribute namespace. So we have to detect and
    // handle this case.
    if (XMLString::equals(qn, XMLUni::fgXMLNSString))
       ns = XMLUni::fgXMLNSURIName;

    cur_->setAttributeNS(ns, qn, attr.getValue (i));
    //std::cerr << XMLString::transcode(attr.getValue(i)) << std::endl;
  }
  depth_++;
}

void parser_impl::endElement (const XMLCh* const /*uri*/,
    const XMLCh* const /*lname*/, const XMLCh* const /*qname*/) {
  // We have an element parent only on depth 2 or greater.

  //  doc_.reset ();
  --depth_;

  if (doc_.get () != 0 && depth_ >= docDepth_ )
    cur_ = static_cast<DOMElement*> (cur_->getParentNode());
}

#if _XERCES_VERSION >= 30000
void parser_impl::characters (const XMLCh* const s, const XMLSize_t length) {
  const XMLCh empty[] = {chNull};

  // Ignore text content (presumably whitespaces) in the root element.
  if (depth_ > 1) {
    DOMText* t = doc_->createTextNode(empty);
    static_cast<DOMTextImpl*> (t)->appendData(s, length);
    cur_->appendChild(t);
  }
}
#else
void parser_impl::characters (const XMLCh* const s, const unsigned int length) {
  // Ignore text content (presumably whitespaces) in the root element.
  if (depth_ > 1) {
    // For Xerces-C++ 2-series we have to make copy.
    xsd::cxx::auto_array<XMLCh> tmp(new XMLCh[length + 1]);
    XMLString::copyNString(tmp.get(), s, length);
    cur_->appendChild(doc_->createTextNode(tmp.get()));
  }
}
#endif

//
// parser
//

parser::parser() : impl_ (new parser_impl) {}
parser::~parser() {}

xml::dom::auto_ptr<DOMDocument> parser::start(istream& is, const string& id, 
    bool val, string schemaDefinition, string schema_major, string schema_minor,
    string schemaNamespace, bool noNameSpace) {
  return impl_->start(is, id, val, schemaDefinition, schema_major, schema_minor, 
                      schemaNamespace, noNameSpace);
}

xml::dom::auto_ptr<DOMDocument> parser::next() {
  return impl_->next ();
}
