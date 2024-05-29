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


void FragSpectrumScanDatabaseBoostdb::terminate()
{
  
  if(bdb)
  {
    bdb->clear();
    delete bdb;
  }
  bdb = 0;
}

std::unique_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseBoostdb::getFSS( unsigned int scanNr ) 
{
  mapdb::const_iterator it;
  it = bdb->find(scanNr);
  if(it == bdb->end()){
    return std::unique_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
  }
  std::istringstream istr (it->second);
  binary_iarchive ia (istr);
  xml_schema::istream<binary_iarchive> is (ia);
  std::unique_ptr< ::percolatorInNs::fragSpectrumScan> ret (new ::percolatorInNs::fragSpectrumScan (is)); 
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
    std::unique_ptr< ::percolatorInNs::fragSpectrumScan> fss (new ::percolatorInNs::fragSpectrumScan (is));
    ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
  }

}

void FragSpectrumScanDatabaseBoostdb::printTab(ostream &tabOutputStream) {
  mapdb::const_iterator it;
  for (it = bdb->begin(); it != bdb->end(); it++) {
    std::istringstream istr (it->second);
    binary_iarchive ia (istr);
    xml_schema::istream<binary_iarchive> is (ia);
    std::unique_ptr< ::percolatorInNs::fragSpectrumScan> fss (new ::percolatorInNs::fragSpectrumScan (is));
    printTabFss(fss, tabOutputStream);
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
