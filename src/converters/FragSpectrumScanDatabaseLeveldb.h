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

#ifndef FRAGSPECTRUMSCANDATABASELEVELDB_H
#define FRAGSPECTRUMSCANDATABASELEVELDB_H

#include "FragSpectrumScanDatabase.h"
#include "leveldb/db.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

class FragSpectrumScanDatabaseLeveldb : public FragSpectrumScanDatabase
{

public:
    
  FragSpectrumScanDatabaseLeveldb(std::string id = 0);

  virtual ~FragSpectrumScanDatabaseLeveldb();

  virtual std::string toString();
    
  virtual bool init(std::string fileName);
  
  virtual void terminate();
  
  virtual std::unique_ptr< ::percolatorInNs::fragSpectrumScan> deserializeFSSfromBinary( char * value, int valueSize );
  
  virtual std::unique_ptr< ::percolatorInNs::fragSpectrumScan> getFSS( unsigned int scanNr );
  
  virtual void print(serializer & ser);
  
  virtual void printTab(ostream &tabOutputStream);
  
  virtual void putFSS( ::percolatorInNs::fragSpectrumScan & fss );
  
private:
  
  XDR xdr;
  xml_schema::buffer buf;
  std::unique_ptr< xml_schema::ostream<XDR> > oxdrp;
  leveldb::DB* bdb;
  leveldb::Options options;
  
};

#endif // FRAGSPECTRUMSCANDATABASELEVELDB_H
