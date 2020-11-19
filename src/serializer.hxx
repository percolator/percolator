// file      : examples/cxx/tree/streaming/serializer.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#ifndef SERIALIZER_HXX
#define SERIALIZER_HXX

#include <string>
#include <iosfwd>
#include <memory> // std::auto_ptr

#include <xercesc/dom/DOMElement.hpp>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

class serializer_impl;

class serializer
{
public:
  ~serializer ();
  serializer ();

  // Start the serialization process.
  //
  void
  start (std::ostream& is, const std::string& encoding = "UTF-8");

  // Serialize next object model fragment into an element with the specified
  // name.
  //
  template <typename T>
  void
  next (const std::string& name, const T& x);

  // Serialize next object model fragment into an element with the specified
  // namespace and qualified name.
  //
  template <typename T>
  void
  next (const std::string& ns, const std::string& name, const T& x);

private:
  serializer (const serializer&);

  serializer&
  operator= (const serializer&);

private:
  xercesc::DOMElement*
  create (const std::string& name);

  xercesc::DOMElement*
  create (const std::string& ns, const std::string& name);

  void
  serialize (xercesc::DOMElement&);

private:
  std::auto_ptr<serializer_impl> impl_;
};

template <typename T>
inline void serializer::
next (const std::string& name, const T& x)
{
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMElement> e (create (name));
  *e << x;
  serialize (*e);
}

template <typename T>
inline void serializer::
next (const std::string& ns, const std::string& name, const T& x)
{
  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMElement> e (create (ns, name));
  *e << x;
  serialize (*e);
}

#endif // SERIALIZER_HXX
