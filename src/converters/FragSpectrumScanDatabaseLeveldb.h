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

  FragSpectrumScanDatabaseLeveldb(const FragSpectrumScanDatabase& original);

  FragSpectrumScanDatabaseLeveldb();

  FragSpectrumScanDatabaseLeveldb(const FragSpectrumScanDatabaseLeveldb& other);

  virtual ~FragSpectrumScanDatabaseLeveldb();

  virtual FragSpectrumScanDatabaseLeveldb& operator=(const FragSpectrumScanDatabaseLeveldb& other);

  virtual bool operator==(const FragSpectrumScanDatabaseLeveldb& other) const;
  
  virtual bool init(std::string fileName);
  
  virtual void terminte();
  
  virtual std::auto_ptr< ::percolatorInNs::fragSpectrumScan> deserializeFSSfromBinary( char * value, int valueSize );
  
  virtual std::auto_ptr< ::percolatorInNs::fragSpectrumScan> getFSS( unsigned int scanNr );
  
  virtual void print(serializer & ser);
  
  virtual void putFSS( ::percolatorInNs::fragSpectrumScan & fss );
  
private:
  
  XDR xdr;
  xml_schema::buffer buf;
  std::auto_ptr< xml_schema::ostream<XDR> > oxdrp;
  leveldb::DB* bdb;
  leveldb::Options options;
  
  
};

#endif // FRAGSPECTRUMSCANDATABASELEVELDB_H
