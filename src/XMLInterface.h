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
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <memory>
#include <cassert>
#include <boost/assign.hpp>

#include "Globals.h"
#include "SetHandler.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "Scores.h"
#include "ProteinProbEstimator.h"
#include "SanityCheck.h"

#ifdef XML_SUPPORT
  #include "Enzyme.h"
  #include "MassHandler.h"
  #include "PseudoRandom.h"
  #include "FeatureMemoryPool.h"
  #include "PSMDescription.h"
  #include "PSMDescriptionDOC.h"
  
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
  XMLInterface(const std::string& xmlOutputFN, const std::string& xmlPeptideOutputFN, const bool xmlSchemaValidation,
               bool printDecoys, bool printExpMass);
  ~XMLInterface();
  
  inline void setXmlOutputFN(std::string outputFN) { xmlOutputFN_ = outputFN; }
  inline std::string getXmlOutputFN() { return xmlOutputFN_; }

  inline void setxmlPeptideOutputFN(std::string outputFN) { xmlPeptideOutputFN_ = outputFN; }
  inline std::string getxmlPeptideOutputFN() { return xmlPeptideOutputFN_; }
  
  inline void setSchemaValidation(bool on) { schemaValidation_ = on; }
  inline void setPrintDecoys(bool decoysOut) { 
    printDecoys_ = decoysOut; 
  }
  inline bool getPrintDecoys() { return printDecoys_; }
  inline void setPrintExpMass(bool printExpMass) { 
    printExpMass_ = printExpMass; 
  }
  inline bool getPrintExpMass() { return printExpMass_; }
  
  int readPin(istream& dataStream, const std::string& xmlInputFN, 
    SetHandler& setHandler, SanityCheck*& pCheck, 
    ProteinProbEstimator* protEstimator, Enzyme*& enzyme);
  int readAndScorePin(istream& dataStream, std::vector<double>& rawWeights, 
    Scores& allScores, const std::string& xmlInputFN,
    SetHandler& setHandler, SanityCheck*& pCheck, 
    ProteinProbEstimator* protEstimator, Enzyme*& enzyme);
  
  void writeXML_PSMs(Scores& fullset);
  void writePeptideXML_PSMs(Scores& fullset, double selectionFdr_);
  void writeXML_Peptides(Scores& fullset);
  void writeXML_Proteins(ProteinProbEstimator* protEstimator);
  void writeXML(Scores& fullset, ProteinProbEstimator* protEstimator, 
                std::string call);
  void writePeptideXML(Scores& fullset, ProteinProbEstimator* protEstimator, 
                std::string call);
 protected:
  map<char, float> getRoughAminoWeightDict();
  
  std::string xmlOutputFN_; 
  std::string xmlPeptideOutputFN_;
  bool schemaValidation_;
  std::string otherCall_;
  
  bool printDecoys_, printExpMass_;
  
  std::string xmlOutputFN_PSMs;
  std::string xmlpeptideOutputFN_PSMs;
  std::string xmlOutputFN_Peptides;
  std::string xmlOutputFN_Proteins;
  
  bool reportUniquePeptides_;
  bool reportPeptideXML_;
  double pi0Psms_;
  double pi0Peptides_;
  unsigned int numberQpsms_;
  
#ifdef XML_SUPPORT
  PSMDescription* readPsm(const ::percolatorInNs::peptideSpectrumMatch &psm, 
                          unsigned scanNumber, bool readProteins,
                          FeatureMemoryPool& featurePool);
  ScanId getScanId(const percolatorInNs::peptideSpectrumMatch& psm, 
                   unsigned scanNumber);
  std::string decoratePeptide(const ::percolatorInNs::peptideType& peptide);
#endif //XML_SUPPORT
  
};

#endif /*XMLINTERFACE_H_*/
