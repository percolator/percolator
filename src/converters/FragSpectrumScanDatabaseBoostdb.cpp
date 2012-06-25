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


#include "FragSpectrumScanDatabaseBoostdb.h"

FragSpectrumScanDatabaseBoostdb::FragSpectrumScanDatabaseBoostdb(std::string id):FragSpectrumScanDatabase(id)
{
  bdb = new mapdb();
}

FragSpectrumScanDatabaseBoostdb::~FragSpectrumScanDatabaseBoostdb()
{

}

string FragSpectrumScanDatabaseBoostdb::toString()
{
  return std::string("FragSpectrumScanDatabaseBoostdb");
}

bool FragSpectrumScanDatabaseBoostdb::init(std::string fileName) {

  //TODO should check the hash table is fine
  bool ret = true;
  return ret;
}


void FragSpectrumScanDatabaseBoostdb::terminte()
{
  bdb->clear();
  delete bdb;
}

std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseBoostdb::getFSS( unsigned int scanNr ) 
{
  mapdb::const_iterator it;
  it = bdb->find(scanNr);
  if(it == bdb->end()){
    return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
  }
  std::istringstream istr (it->second);
  binary_iarchive ia (istr);
  xml_schema::istream<binary_iarchive> is (ia);
  std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret (new ::percolatorInNs::fragSpectrumScan (is)); 
  return ret;      
}

void FragSpectrumScanDatabaseBoostdb::print(serializer & ser) 
{

  mapdb::const_iterator it;
  for (it = bdb->begin(); it != bdb->end(); it++) 
  {
    std::istringstream istr (it->second);
    binary_iarchive ia (istr);
    xml_schema::istream<binary_iarchive> is (ia);
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss (new ::percolatorInNs::fragSpectrumScan (is));
    ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
  }

}

void FragSpectrumScanDatabaseBoostdb::putFSS( ::percolatorInNs::fragSpectrumScan & fss ) 
{
  std::ostringstream ostr;
  binary_oarchive oa (ostr);
  xml_schema::ostream<binary_oarchive> os (oa);
  os << fss;
  ostr.flush();
  ::percolatorInNs::fragSpectrumScan::scanNumber_type key = fss.scanNumber();
  (*bdb)[key] = ostr.str();
  ostr.str(""); // reset the string
}
