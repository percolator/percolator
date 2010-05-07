// file      : examples/cxx/tree/streaming/parser.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#ifndef PARSER_HXX
#define PARSER_HXX

#include <string>
#include <iosfwd>
#include <memory> // std::auto_ptr

#include <xercesc/dom/DOMDocument.hpp>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

class parser_impl;

class parser
{
public:
  ~parser ();
  parser ();

  // The start function returns a "carcase" of the complete document. That
  // is, the root element with all the attributes but without any content.
  //
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument>
  start (std::istream& is, const std::string& id, bool validate);

  // The next function returns next first-level element with all its
  // attributes and content or 0 if no more available.
  //
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument>
  next ();

private:
  parser (const parser&);

  parser&
  operator= (const parser&);

private:
  std::auto_ptr<parser_impl> impl_;
};

#endif // PARSER_HXX
