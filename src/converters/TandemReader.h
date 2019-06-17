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


#ifndef TANDEMREADER_H
#define TANDEMREADER_H

#include "Reader.h"
#include "parser.hxx"
#include "tandem2011.12.01.1.hxx"
#include "FragSpectrumScanDatabase.h"
#include <boost/foreach.hpp>

using namespace std;
using namespace xercesc;

typedef map<std::string, set<std::string> > peptideProteinMapType;

class TandemReader: public Reader {

 public:
  
  TandemReader(ParseOptions po);
  
  virtual ~TandemReader();
  
  void read(const std::string &fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database);
  
  bool checkValidity(const std::string &file);
  
  bool checkIsMeta(const std::string &file);
 
  void getMaxMinCharge(const std::string &fn, bool isDecoy);
  
  void addFeatureDescriptions(bool doEnzyme);
  
 protected:
  
  //Variables
  bool x_score, y_score, z_score, a_score, b_score, c_score;
  bool firstPSM;
  
  //Functions
  void readSpectra(const tandem_ns::group &groupObj,bool isDecoy,
		  boost::shared_ptr<FragSpectrumScanDatabase> database,const std::string &fn);
  
  void getPeptideProteinMap(const tandem_ns::group &groupObj,
      peptideProteinMapType &peptideProteinMap);
  
  void createPSM(const tandem_ns::peptide::domain_type &domain,double parenIonMass,unsigned charge,
		  double sumI,double maxI,bool isDecoy, boost::shared_ptr<FragSpectrumScanDatabase> database,
		  const peptideProteinMapType &peptideProteinMap,const string &psmId, int spectraId);
  
  static const std::map<string,double> tandemFeaturesDefaultValue;
};

#endif //TANDEMREADER_H
