#include "FragSpectrumScanDatabaseTokyodb.h"

extern "C"
typedef  void (*xdrrec_create_p) (
    XDR*,
    unsigned int write_size,
    unsigned int read_size,
    void* user_data,
    int (*read) (void* user_data, char* buf, int n),
    int (*write) (void* user_data, char* buf, int n));



FragSpectrumScanDatabaseTokyoDB::FragSpectrumScanDatabaseTokyoDB(std::string id):FragSpectrumScanDatabase(id)
{
  bdb = 0;
  xdrrec_create_p xdrrec_create_ = reinterpret_cast<xdrrec_create_p> (::xdrrec_create);
  xdrrec_create_ (&xdr, 0, 0, reinterpret_cast<char*> (&buf), 0, &overflow);
  xdr.x_op = XDR_ENCODE;
  std::auto_ptr< xml_schema::ostream<XDR> > tmpPtr(new xml_schema::ostream<XDR>(xdr)) ;
  assert(tmpPtr.get());
  oxdrp=tmpPtr;
}

FragSpectrumScanDatabaseTokyoDB::~FragSpectrumScanDatabaseTokyoDB()
{

}

std::string FragSpectrumScanDatabaseTokyoDB::toString()
{
  return std::string("FragSpectrumScanDatabaseTokyoDB");
}

//   The function "tcbdbopen" in Tokyo Cabinet does not have O_EXCL as is
//   possible in the unix system call open (see "man 2 open"). This may be a
//   security issue if the filename to the Tokyo cabinet database is in a
//   directory that other users have write access to. They could add a symbolic
//   link pointing somewhere else. It would be better if Tokyo Cabinet would
//   fail if the database existed in our case when we use a temporary file.
bool FragSpectrumScanDatabaseTokyoDB::init(std::string fileName) 
{
  bdb = tcbdbnew();
  assert(bdb);
  bool ret =  tcbdbsetcmpfunc(bdb, tccmpint32, NULL);
  assert(ret);
  if(!tcbdbopen(bdb, fileName.c_str(), BDBOWRITER | BDBOTRUNC | BDBOREADER | BDBOCREAT )){
    int errorcode = tcbdbecode(bdb);
    fprintf(stderr, "open error: %s\n", tcbdberrmsg(errorcode));
    exit(EXIT_FAILURE);
  }
  ret = unlink( fileName.c_str() );
  assert(! ret);
  return ret;
}

void FragSpectrumScanDatabaseTokyoDB::terminte()
{
  tcbdbdel(bdb);   
}


std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseTokyoDB::deserializeFSSfromBinary( char * value, int valueSize ) 
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
  std::auto_ptr< percolatorInNs::fragSpectrumScan> fss (new percolatorInNs::fragSpectrumScan(ixdr));
  
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

std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabaseTokyoDB::getFSS( unsigned int scanNr ) 
{
  assert(bdb);
  int valueSize = 0;
  char * value = ( char * ) tcbdbget(bdb, ( const char* ) &scanNr, sizeof( scanNr ), &valueSize);
  if(!value) {
    return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
  }
  std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret(deserializeFSSfromBinary(value,valueSize));  
  free(value);
  return ret;
}

void FragSpectrumScanDatabaseTokyoDB::print(serializer & ser) 
{
  BDBCUR *cursor;
  char *key;
  assert(bdb);
  cursor = tcbdbcurnew(bdb);
  assert(cursor);
  tcbdbcurfirst(cursor);
  // using tcbdbcurkey3 is probably faster
  int keySize;
  int valueSize;
  while (( key = static_cast< char * > ( tcbdbcurkey(cursor,&keySize)) ) != 0 ) 
  {
    char * value = static_cast< char * > ( tcbdbcurval(cursor,&valueSize));
    if(value)
    {
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss(deserializeFSSfromBinary(value,valueSize));
      ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
      free(value);
    }
    free(key);
    tcbdbcurnext(cursor);
  }
  tcbdbcurdel(cursor);
}

void FragSpectrumScanDatabaseTokyoDB::putFSS( ::percolatorInNs::fragSpectrumScan & fss ) 
{   
  assert(bdb);
  ::percolatorInNs::fragSpectrumScan::scanNumber_type key = fss.scanNumber();
  *oxdrp << fss;
  xdrrec_endofrecord (&xdr, true);
  size_t keySize = sizeof(key);
  size_t valueSize(buf.size ());
  if(!tcbdbput(bdb, ( const char * ) &key, keySize, buf.data (), buf.size () ))
  {
    int  errorcode = tcbdbecode(bdb);
    fprintf(stderr, "put error: %s\n", tcbdberrmsg(errorcode));
    exit(EXIT_FAILURE);
  }
  buf.size(0);
}
