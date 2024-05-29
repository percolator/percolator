#include "FragSpectrumScanDatabaseLeveldb.h"


extern "C"
typedef  void (*xdrrec_create_p) (
    XDR*,
    unsigned int write_size,
    unsigned int read_size,
    void* user_data,
    int (*read) (void* user_data, char* buf, int n),
    int (*write) (void* user_data, char* buf, int n));


FragSpectrumScanDatabaseLeveldb::FragSpectrumScanDatabaseLeveldb(std::string id):FragSpectrumScanDatabase(id)
{
  bdb = 0;
  xdrrec_create_p xdrrec_create_ = reinterpret_cast<xdrrec_create_p> (::xdrrec_create);
  xdrrec_create_ (&xdr, 0, 0, reinterpret_cast<char*> (&buf), 0, &overflow);
  xdr.x_op = XDR_ENCODE;
  std::unique_ptr< xml_schema::ostream<XDR> > tmpPtr(new xml_schema::ostream<XDR>(xdr)) ;
  assert(tmpPtr.get());
  oxdrp=tmpPtr;
}

FragSpectrumScanDatabaseLeveldb::~FragSpectrumScanDatabaseLeveldb()
{

}

std::string FragSpectrumScanDatabaseLeveldb::toString()
{
  return std::string("FragSpectrumScanDatabaseLeveldb");
}

bool FragSpectrumScanDatabaseLeveldb::init(std::string fileName) 
{
  options.create_if_missing = true;
  options.error_if_exists = true;
  options.max_open_files = 100;
  options.write_buffer_size = 4194304*2; //8 MB
  options.block_size = 4096*4; //16K
  leveldb::Status status = leveldb::DB::Open(options, fileName.c_str(), &bdb);
  if (!status.ok()){ 
    std::cerr << status.ToString() << endl;
  }
  bool ret = status.ok();
  return ret;
}

void FragSpectrumScanDatabaseLeveldb::terminate()
{
  if(bdb) delete(bdb);
  bdb = 0;
}

std::unique_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseLeveldb::deserializeFSSfromBinary( char * value, int valueSize ) 
{
  xml_schema::buffer buf2;
  buf2.capacity(valueSize);
  memcpy(buf2.data(), value, valueSize);
  buf2.size(valueSize);
  underflow_info ui;
  ui.buf = &buf2;
  ui.pos = 0;
  XDR xdr2;
  xdrrec_create_p xdrrec_create_ = reinterpret_cast<xdrrec_create_p> (::xdrrec_create);
  xdrrec_create_ (&xdr2, 0, 0, reinterpret_cast<char*> (&ui), &underflow, 0);
  xdr2.x_op = XDR_DECODE;
  xml_schema::istream<XDR> ixdr(xdr2);
  xdrrec_skiprecord(&xdr2);
  std::unique_ptr< percolatorInNs::fragSpectrumScan> fss (new percolatorInNs::fragSpectrumScan(ixdr));
  
  //TODO this gives too many arguments in MINGW
  //xdr_destroy (&xdr2);
  if((&xdr2)->x_ops->x_destroy)			
  {
#if defined __MINGW__ or defined __WIN32__
    (*(&xdr2)->x_ops->x_destroy);
#else
    (*(&xdr2)->x_ops->x_destroy)(&xdr2);
#endif
  }
  return fss;
}

std::unique_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseLeveldb::getFSS( unsigned int scanNr ) 
{
  assert(bdb);
  std::string skey = boost::lexical_cast<std::string>(scanNr);
  leveldb::Slice s1(skey);
  leveldb::Iterator* itr = bdb->NewIterator(leveldb::ReadOptions());
  itr->Seek(s1);
  if(!itr->Valid() || s1 != itr->key()){
    delete itr;
    return std::unique_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
  }
  char *retvalue = const_cast<char*>(itr->value().data());
  std::unique_ptr< ::percolatorInNs::fragSpectrumScan> ret(deserializeFSSfromBinary(retvalue,itr->value().size()));
  delete itr;
  return ret;
}

void FragSpectrumScanDatabaseLeveldb::print(serializer & ser) 
{
  assert(bdb);
  leveldb::Iterator* it = bdb->NewIterator(leveldb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    char *retvalue = const_cast<char*>(it->value().data());
    std::unique_ptr< ::percolatorInNs::fragSpectrumScan> fss(deserializeFSSfromBinary(retvalue,it->value().size()));
    ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
  }
  delete it;

}

void FragSpectrumScanDatabaseLeveldb::printTab(ostream &tabOutputStream) {
  assert(bdb);
  leveldb::Iterator* it = bdb->NewIterator(leveldb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    char *retvalue = const_cast<char*>(it->value().data());
    std::unique_ptr< ::percolatorInNs::fragSpectrumScan> fss(deserializeFSSfromBinary(retvalue,it->value().size()));
    printTabFss(fss, tabOutputStream);
  }
  delete it;
}


void FragSpectrumScanDatabaseLeveldb::putFSS( ::percolatorInNs::fragSpectrumScan & fss ) 
{
  assert(bdb);
  *oxdrp << fss;
  xdrrec_endofrecord (&xdr, true);
  leveldb::WriteOptions write_options;
  ::percolatorInNs::fragSpectrumScan::scanNumber_type key = fss.scanNumber();
  std::string skey = boost::lexical_cast<std::string>(key);
  leveldb::Slice s2(buf.data(),buf.size());
  leveldb::Slice s1(skey);
  leveldb::Status status = bdb->Put(write_options,s1,s2);
  if(!status.ok())
  {
    throw MyException(status.ToString());
  }
  buf.size(0);
}
