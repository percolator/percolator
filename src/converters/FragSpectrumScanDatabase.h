#ifndef FRAGSPECTRUMSCANDATABASE_H
#define FRAGSPECTRUMSCANDATABASE_H

#include <memory>   // std::auto_ptr
#include <tcbdb.h>
#include <rpc/xdr.h>
#include <rpc/types.h>
#include <iostream>
#include <xercesc/dom/DOM.hpp>
#include <cstddef>  // size_t
#include <cstring>  // memcpy
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <map>
#include <string>
#include <memory>    // std::auto_ptr
#include <iostream>
#include "config.h"
#include "serializer.hxx"
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include <xercesc/util/XMLUni.hpp>
#include <iostream>
#include <xercesc/dom/DOM.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>
#include <xsd/cxx/xml/dom/serialization-source.hxx>
#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xercesc/dom/DOMElement.hpp>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include "percolator_in.hxx"

using namespace std;
using namespace percolatorInNs;

class serializer;

class FragSpectrumScanDatabase {
  public:
    FragSpectrumScanDatabase(string id=0);
    ~FragSpectrumScanDatabase(){}
    FragSpectrumScanDatabase(const FragSpectrumScanDatabase& original){}
    bool init(std::string filename);
    bool initRTime(map<int, vector<double> >* scan2rt_par);
    auto_ptr<fragSpectrumScan> getFSS( unsigned int scanNr );
    auto_ptr<fragSpectrumScan> deserializeFSSfromBinary(char* value,
        int valueSize);
    void putFSS(fragSpectrumScan & fss );
    void savePsm(unsigned int scanNr, auto_ptr<peptideSpectrumMatch> psm_p );
    void print(serializer & ser );
    void terminte();
    string id;
  protected:
    XDR xdr;
    xml_schema::buffer buf;
    TCBDB* bdb;
    // pointer to retention times
    map<int, vector<double> >* scan2rt;
    // is scoped_ptr possible here?
    std::auto_ptr< xml_schema::ostream<XDR> > oxdrp;
};

#endif
