// file      : examples/cxx/tree/streaming/parser.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#ifndef PARSER_HXX
#define PARSER_HXX

#include <string>
#include <fstream>
#include <iosfwd>
#include <memory> // std::auto_ptr
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include "Globals.h"
using namespace xercesc;

static const XMLCh spectrumIdentificationListStr[] = { 
  chLatin_S, chLatin_p, chLatin_e, chLatin_c, chLatin_t,chLatin_r, chLatin_u, chLatin_m, 
  chLatin_I, chLatin_d, chLatin_e, chLatin_n, chLatin_t, chLatin_i, chLatin_f, chLatin_i, chLatin_c, chLatin_a, chLatin_t, chLatin_i, chLatin_o, chLatin_n,
  chLatin_L, chLatin_i, chLatin_s, chLatin_t, chNull };

static const XMLCh spectrumIdentificationResultStr[] = { 
  chLatin_S, chLatin_p, chLatin_e, chLatin_c, chLatin_t, chLatin_r, chLatin_u, chLatin_m, 
  chLatin_I, chLatin_d, chLatin_e, chLatin_n, chLatin_t, chLatin_i, chLatin_f, chLatin_i, chLatin_c, chLatin_a, chLatin_t, chLatin_i, chLatin_o, chLatin_n, 
  chLatin_R, chLatin_e, chLatin_s, chLatin_u, chLatin_l, chLatin_t, chNull };

class parser_impl;

class parser {
 public:
  ~parser();
  parser();

  // The start function returns a "carcase" of the complete document. That
  // is, the root element with all the attributes but without any content.
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> start(std::istream& is, 
    const std::string& id, bool validate, std::string schemaDefinition, 
    std::string schema_major, std::string schema_minor,
    std::string schemaNamespace = "http://per-colator.com/percolator_in/", 
    bool noNameSpace = false);

  // The next function returns next first-level element with all its
  // attributes and content or 0 if no more available.
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> next();

 private:
  parser(const parser&);

  parser& operator= (const parser&);

 private:
  std::auto_ptr<parser_impl> impl_;
};

#endif // PARSER_HXX
