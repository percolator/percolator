/*
    Copyright 2012 <copyright holder> <email>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#ifndef MSGFDBREADER_H
#define MSGFDBREADER_H

#include "Reader.h"

class msgfdbReader: public Reader
{

public:
  
  msgfdbReader();
  msgfdbReader(const msgfdbReader& other);
  virtual ~msgfdbReader();
  virtual msgfdbReader& operator=(const msgfdbReader& other);
  virtual bool operator==(const msgfdbReader& other) const;

  void read(const std::string fn,
    ::percolatorInNs::featureDescriptions& fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool is_decoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType t, FragSpectrumScanDatabase* database);
  
private:
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);
  std::string remove_endl(std::string s );
  
  void readPSM(std::string line,bool isDecoy,
    const ParseOptions& po,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss,
     int minCharge, int maxCharge, std::string fileId,
     FragSpectrumScanDatabase* database);
  
  void searchMaxMinCharge(const std::string fn, int* maxCharge,int* minCharge);
  
  double calculatePepMAss(std::string pepsequence,double charge);

};

#endif //MSGFDBREADER_H
