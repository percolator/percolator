/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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

#ifndef MZIDENTMLREADER_H
#define MZIDENTMLREADER_H

#include <Reader.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <map>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <boost/foreach.hpp>
#include <algorithm>
#include "Enzyme.h"
#include "DataSet.h"
#include "MassHandler.h"
#include "FragSpectrumScanDatabase.h"
#include "parser.hxx"
#include "mzIdentML1.0.0.hxx"

using namespace std;
using namespace xercesc;
typedef map<std::string, mzIdentML_ns::SequenceCollectionType::Peptide_type *> peptideMapType;
typedef map<std::string, int> scanNumberMapType;

class MzidentmlReader: public Reader
{

public:

  MzidentmlReader(ParseOptions po);
  virtual ~MzidentmlReader();

  virtual void read(const std::string fn, percolatorInNs::featureDescriptions& fds, 
		    percolatorInNs::experiment::fragSpectrumScan_sequence& fsss, bool is_decoy, 
		    FragSpectrumScanDatabase* database);
  
  virtual bool checkValidity(std::string file);
 
  virtual void getMaxMinCharge(std::string fn);
  
  int loadFromTargetOrDecoyFile( const char * fileName,  bool isDecoy,  
			       percolatorInNs::featureDescriptions & fdesFirstFile, 
			       FragSpectrumScanDatabase & database,  int * scanNumber, 
			       scanNumberMapType & scanNumberMap  );
 
  void createPSM(const ::mzIdentML_ns::PSI_PI_analysis_search_SpectrumIdentificationItemType & item, 
		peptideMapType & peptideMap,::percolatorInNs::fragSpectrumScan::experimentalMassToCharge_type 
		experimentalMassToCharge,bool isDecoy, percolatorInNs::featureDescriptions & fdesFirstFile,  
		::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_sequence & psm_sequence );
  
};

#endif // MZIDENTMLREADER_H
