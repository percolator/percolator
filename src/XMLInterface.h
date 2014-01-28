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
#ifndef XMLINTERFACE_H_
#define XMLINTERFACE_H_

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Globals.h"
#include "SetHandler.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "Scores.h"
#include "ProteinProbEstimator.h"
#ifdef XML_SUPPORT
  #include "SanityCheck.h"
  #include "SqtSanityCheck.h"
  #include "FidoInterface.h"
  #include "Enzyme.h"
  
  #include "parser.hxx"
  #include "serializer.hxx"
  #include <xercesc/dom/DOM.hpp>
  #include <xercesc/util/XMLString.hpp>
  #include <xercesc/parsers/XercesDOMParser.hpp>
  #include <xercesc/sax/HandlerBase.hpp>
  #include <xercesc/util/PlatformUtils.hpp>
  #include "percolator_in.hxx"
#endif //XML_SUPPORT

class XMLInterface {
  
  public:
    XMLInterface();
    
    inline void setXmlInputFN(std::string inputFN) { xmlInputFN = inputFN; }
    inline std::string getXmlInputFN() { return xmlInputFN; }
    inline void setXmlOutputFN(std::string outputFN) { xmlOutputFN = outputFN; }
    inline std::string getXmlOutputFN() { return xmlOutputFN; }
    
    inline void setSchemaValidation(bool on) { schemaValidation = on; }
    
    int readPin(SetHandler & setHandler, SanityCheck *& pCheck, ProteinProbEstimator * protEstimator);
    
    inline void setPi0Peptides(double pi0) { pi_0_peptides = pi0; }
    inline void setPi0Psms(double pi0) { pi_0_psms = pi0; }
    inline void setNumberQpsms(double nq) { numberQpsms = nq; }
    
    void writeXML_PSMs(Scores & fullset);
    void writeXML_Peptides(Scores & fullset);
    void writeXML_Proteins(ProteinProbEstimator * protEstimator);
    void writeXML(Scores & fullset, ProteinProbEstimator * protEstimator, std::string call);
    
  protected:
    std::string xmlInputFN;
    bool schemaValidation;
    std::string otherCall;
    
    std::string xmlOutputFN;
    std::string xmlOutputFN_PSMs;
    std::string xmlOutputFN_Peptides;
    std::string xmlOutputFN_Proteins;
    
    bool reportUniquePeptides;
    double pi_0_psms;
    double pi_0_peptides;
    double numberQpsms;
    
};

#ifdef XML_SUPPORT

/** some constant strings to be used to compare xml strings **/

//databases
static const XMLCh databasesStr[] = {
  xercesc::chLatin_d, xercesc::chLatin_a, xercesc::chLatin_t, xercesc::chLatin_a, 
  xercesc::chLatin_b, xercesc::chLatin_a, xercesc::chLatin_s, xercesc::chLatin_e, 
  xercesc::chLatin_s, xercesc::chNull };
      
//calibration
static const XMLCh calibrationStr[] = { 
  xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_l, xercesc::chLatin_i, 
  xercesc::chLatin_b, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_t, 
  xercesc::chLatin_i, xercesc::chLatin_o, xercesc::chLatin_n, xercesc::chNull };

//proteins
static const XMLCh proteinsStr[] = { 
  xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t, 
  xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chLatin_s, 
  xercesc::chNull };

//protein      
static const XMLCh proteinStr[] = { 
  xercesc::chLatin_p, xercesc::chLatin_r, xercesc::chLatin_o, xercesc::chLatin_t, 
  xercesc::chLatin_e, xercesc::chLatin_i, xercesc::chLatin_n, xercesc::chNull };
      
//fragSpectrumScan 
static const XMLCh fragSpectrumScanStr[] = { 
  xercesc::chLatin_f, xercesc::chLatin_r, xercesc::chLatin_a, xercesc::chLatin_g, 
  xercesc::chLatin_S, xercesc::chLatin_p, xercesc::chLatin_e, xercesc::chLatin_c, 
  xercesc::chLatin_t, xercesc::chLatin_r, xercesc::chLatin_u, xercesc::chLatin_m,
  xercesc::chLatin_S, xercesc::chLatin_c, xercesc::chLatin_a, xercesc::chLatin_n, 
  xercesc::chNull };  
      
#endif //XML_SUPPORT


#endif /*XMLINTERFACE_H_*/
