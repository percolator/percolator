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

#ifndef MZIDENTMLREADER_H
#define MZIDENTMLREADER_H

#include <Reader.h>
#include "FragSpectrumScanDatabase.h"
#include "parser.hxx"
#include "mzIdentML1.1.0.hxx"

using namespace std;
using namespace xercesc;
typedef map<std::string, mzIdentML_ns::SequenceCollectionType::Peptide_type *> peptideMapType;
typedef map<std::string, mzIdentML_ns::SequenceCollectionType::DBSequence_type *> proteinMapType;
typedef multimap<std::string, mzIdentML_ns::PeptideEvidenceType *> peptideEvidenceMapType_peptideid;
typedef map<std::string, mzIdentML_ns::PeptideEvidenceType *> peptideEvidenceMapType;
typedef map<std::string, int> scanNumberMapType;

struct RetrieveValue
{
  template <typename T>
  typename T::second_type operator()(T keyValuePair) const
  {
    return keyValuePair.second;
  }
};

class MzidentmlReader: public Reader
{

public:

  MzidentmlReader(ParseOptions po);

  virtual ~MzidentmlReader();

  void read(const std::string &fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database);

  virtual bool checkValidity(const std::string &file) = 0;

  bool checkIsMeta(const std::string &file);

  void getMaxMinCharge(const std::string &fn, bool isDecoy);
  virtual void searchEngineSpecificParsing(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, int itemCount);

  virtual void addFeatureDescriptions(bool doEnzyme) = 0;

  virtual void createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
		  ::percolatorInNs::fragSpectrumScan::experimentalMass_type experimentalMass,
		   bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database,
		   const std::string & fn) = 0;

  void cleanHashMaps();

 protected :

    peptideMapType peptideMap;
    proteinMapType proteinMap;
    peptideEvidenceMapType peptideEvidenceMap;
};

#endif // MZIDENTMLREADER_H
