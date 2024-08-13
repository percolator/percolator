/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef FRAGSPECTRUMSCANDATABASE_H
#define FRAGSPECTRUMSCANDATABASE_H
 
#include <memory>   // std::unique_ptr
#include <iostream>
#include <cstddef>  // size_t
#include <cstring>  // memcpy
#include <map>
#include <list>
#include <string>
#include <algorithm>
#include <cmath>
#include "Globals.h"
#include "MassHandler.h"
#include "serializer.hxx"
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>
#include <xsd/cxx/xml/dom/serialization-source.hxx>
#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include "percolator_in.hxx"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace percolatorInNs;
class serializer;

struct underflow_info
{
    xml_schema::buffer* buf;
    size_t pos;
};

extern "C" int
overflow (void* user_data, char* buf, int n);

extern "C" int
underflow (void* user_data, char* buf, int n);


class FragSpectrumScanDatabase {
  
  public:
    
    FragSpectrumScanDatabase(string id=0);
    
    ~FragSpectrumScanDatabase(){};
    
    bool initRTime(map<int, vector<double> >* scan2rt_par);
    
    void savePsm(unsigned int scanNr, unique_ptr<peptideSpectrumMatch> psm_p );
    
    virtual std::string toString() = 0;
    
    virtual void putFSS(fragSpectrumScan & fss )= 0;
    
    virtual bool init(std::string filename) = 0;
    
    virtual unique_ptr<fragSpectrumScan> getFSS( unsigned int scanNr ) = 0;
    
    virtual unique_ptr<fragSpectrumScan> deserializeFSSfromBinary(char* value,int valueSize) = 0;
    
    virtual void print(serializer & ser ) = 0;
    virtual void printTab(ostream &tabOutputStream) = 0;
    void printTabFss(std::unique_ptr< ::percolatorInNs::fragSpectrumScan>& fss, std::ostream &tabOutputStream); // Pass by reference
    std::string decoratePeptide(const ::percolatorInNs::peptideType& peptide);
    
    virtual void terminate() = 0;
    
    std::string id;
  
  protected:
    // pointer to retention times
    map<int, vector<double> >* scan2rt;
    
    
};

#endif
