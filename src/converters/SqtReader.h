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

#ifndef SQTREADER_H
#define SQTREADER_H

#include "Reader.h"

class SqtReader: public Reader
{

public:
  
  SqtReader(ParseOptions po);
  
  virtual ~SqtReader();
  
  virtual void read(const std::string fn, bool isDecoy,boost::shared_ptr<FragSpectrumScanDatabase> database);

  void readSectionS(std::string record,std::set<int> & theMs, bool isDecoy,
	            std::string psmId,boost::shared_ptr<FragSpectrumScanDatabase> database);

  void readPSM(bool isDecoy, const std::string &in,int match, 
	       std::string psmId,boost::shared_ptr<FragSpectrumScanDatabase> database);
  
  virtual bool checkValidity(const std::string file);
 
  virtual void getMaxMinCharge(std::string fn);
  
  virtual void readRetentionTime(std::string filename);
	
  virtual void storeRetentionTime(boost::shared_ptr<FragSpectrumScanDatabase> database);
  
  virtual void addFeatureDescriptions(bool doEnzyme,const std::string& aaAlphabet,std::string fn);
};

#endif // SQTREADER_H
