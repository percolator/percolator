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


#ifndef SQTREADER_H
#define SQTREADER_H

#include "Reader.h"

class SqtReader: public Reader
{

public:
  
  SqtReader();
  SqtReader(const SqtReader& other);
  virtual ~SqtReader();
  virtual SqtReader& operator=(const SqtReader& other);
  virtual bool operator==(const SqtReader& other) const;

  void read(const std::string fn,
    ::percolatorInNs::featureDescriptions& fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool is_decoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType t, FragSpectrumScanDatabase* database);

  void readSectionS( std::string record,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
     std::set<int> & theMs, bool is_decoy, const ParseOptions& po,
     int minCharge, int maxCharge, std::string psmId,
     FragSpectrumScanDatabase* database);

  void readPSM(bool is_decoy, const std::string &in, int match,
    const ParseOptions& po,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss,
     int minCharge, int maxCharge, std::string psmId,
     FragSpectrumScanDatabase* database );

};

#endif // SQTREADER_H
