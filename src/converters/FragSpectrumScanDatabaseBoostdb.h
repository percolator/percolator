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

#ifndef FRAGSPECTRUMSCANDATABASEBOOSTDB_H
#define FRAGSPECTRUMSCANDATABASEBOOSTDB_H

#if defined(_MSC_VER)
#define _CMATH_IN_CRT 1
#include <float.h>
#endif

#include "FragSpectrumScanDatabase.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using boost::archive::binary_oarchive;
using boost::archive::binary_iarchive;
typedef std::map<unsigned int, std::string, std::less<unsigned int> > mapdb;

class FragSpectrumScanDatabaseBoostdb: public FragSpectrumScanDatabase
{

public:

  FragSpectrumScanDatabaseBoostdb(std::string id = 0);
  
  virtual ~FragSpectrumScanDatabaseBoostdb();
  
  virtual std::string toString();  
  
  virtual bool init(std::string fileName);
  
  virtual void terminate();
  
  virtual std::unique_ptr< ::percolatorInNs::fragSpectrumScan> getFSS( unsigned int scanNr );
  
  virtual void print(serializer & ser);
  
  virtual void printTab(ostream &tabOutputStream);
  
  virtual void putFSS( ::percolatorInNs::fragSpectrumScan & fss );
  
  virtual unique_ptr<fragSpectrumScan> deserializeFSSfromBinary(char* value,int valueSize){throw MyException("deserializeFSSfromBinary is not implemented");};
  
private:
 
   mapdb *bdb;
};

#endif // FRAGSPECTRUMSCANDATABASEBOOSTDB_H
