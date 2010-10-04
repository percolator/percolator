#include <memory>   // std::auto_ptr
#include <iostream>

#include "percolator_in.hxx"

using namespace std;

#include <xercesc/dom/DOM.hpp>
#include <cstddef>  // size_t
#include <cstring>  // memcpy
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <map>
#include <string>
#include <memory>    // std::auto_ptr
#include <iostream>
#include "FragSpectrumScanDatabase.h"
#include "config.h"
#include "serializer.hxx"


#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xercesc/util/XMLUni.hpp>
#include <iostream>
#include <xercesc/dom/DOM.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>
#include <xsd/cxx/xml/dom/serialization-source.hxx>
#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xercesc/dom/DOMElement.hpp>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>


struct underflow_info
{
  xml_schema::buffer* buf;
  size_t pos;
};

extern "C" int
overflow (void* user_data, char* buf, int n);

extern "C" int
underflow (void* user_data, char* buf, int n);



extern "C"
  typedef  void (*xdrrec_create_p) (
				    XDR*,
				    unsigned int write_size,
				    unsigned int read_size,
				    void* user_data,
				    int (*read) (void* user_data, char* buf, int n),
				    int (*write) (void* user_data, char* buf, int n));


/*
extern "C"
typedef  void (*xdrrec_create_p) (
  XDR*,
  unsigned int write_size,
  unsigned int read_size,
  void* user_data,
  int (*read) (void* user_data, char* buf, int n),
  int (*write) (void* user_data, char* buf, int n));
*/

FragSpectrumScanDatabase::FragSpectrumScanDatabase() : bdb(0), scan2rt(0) {
  xdrrec_create_p xdrrec_create_ = reinterpret_cast<xdrrec_create_p> (::xdrrec_create);
   xdrrec_create_ (&xdr, 0, 0, reinterpret_cast<char*> (&buf), 0, &overflow);
   xdr.x_op = XDR_ENCODE;
   std::auto_ptr< xml_schema::ostream<XDR> > tmpPtr(new xml_schema::ostream<XDR>(xdr)) ;
   assert(tmpPtr.get());
   oxdrp=tmpPtr;
}


void FragSpectrumScanDatabase::savePsm( unsigned int scanNr,
    double observedMassCharge,
    std::auto_ptr< percolatorInNs::peptideSpectrumMatch > psm_p ) {

  std::auto_ptr< ::percolatorInNs::fragSpectrumScan>  fss = getFSS( scanNr );
  // if FragSpectrumScan does not yet exist, create it
  if ( ! fss.get() ) {
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fs_p( new ::percolatorInNs::fragSpectrumScan(scanNr, observedMassCharge));
    fss = fs_p;
    // if a retention time has been calculated, include it in the FragSpectrumScan
    if(scan2rt != 0){
      // retrieve retention time
      //double retTime = scan2rt->find(scanNr)->second;
      // fs_p.get()->observedTime().set(retTime);
    }
    //fs_p->observedTime().set(1.0);
  }
  // add the psm to the FragSpectrumScan
  fss->peptideSpectrumMatch().push_back( psm_p );
  putFSS( *fss );
  return;
}

bool FragSpectrumScanDatabase::init(std::string fileName) {
  bdb = tcbdbnew();
  assert(bdb);
/*
  retu = tchdbsetxmsiz(hdb, 200*1024*1024);
  assert( retu );
  retu = tchdbsetcache(hdb, 1000);
  assert( retu );
*/
  bool ret =  tcbdbsetcmpfunc(bdb, tccmpint32, NULL);
  assert(ret);
  if(!tcbdbopen(bdb, fileName.c_str(), BDBOWRITER | BDBOTRUNC | BDBOREADER | BDBOCREAT )){
    int errorcode = tcbdbecode(bdb);
    fprintf(stderr, "open error: %s\n", tcbdberrmsg(errorcode));
    exit(EXIT_FAILURE);
  }
  // unlink => Potential race condition, but Unix seems to lack unlink() with an file descriptor argument.

  ret = unlink( fileName.c_str() );
  assert(! ret);
}

bool FragSpectrumScanDatabase::initRTime(map<int, double>* scan2rt_par) {
  // add pointer to retention times table in sqt2pin (if any)
  scan2rt=scan2rt_par;
}

std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabase::deserializeFSSfromBinary( char * value, int valueSize ) {
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
  xdr_destroy (&xdr2);
  return fss;
}

std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabase::getFSS( unsigned int scanNr ) {
  assert(bdb);
  int valueSize = 0;
  char * value = ( char * ) tcbdbget(bdb, ( const char* ) &scanNr, sizeof( scanNr ), &valueSize);
  if(!value) {
    // ecode = tcbdbecode(bdb);
    // fprintf(stderr, "get error: %s\n", tcbdberrmsg(ecode));
    return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
  }
  std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret(deserializeFSSfromBinary(value,valueSize));
  free(value);
  return ret;
}

void FragSpectrumScanDatabase::print(serializer & ser) {
  BDBCUR *cursor;
  char *key;
  assert(bdb);
  cursor = tcbdbcurnew(bdb);
  assert(cursor);
  tcbdbcurfirst(cursor);
  // using tcbdbcurkey3 is probably faster
  int keySize;
  int valueSize;
  while (( key = static_cast< char * > ( tcbdbcurkey(cursor,&keySize)) ) != 0 ) {
    char * value = static_cast< char * > ( tcbdbcurval(cursor,&valueSize));
    if(value){
       std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss(deserializeFSSfromBinary(value,valueSize));
       // XMLCh * strName = xercesc::XMLString::transcode("XML 1.0 Traversal 2.0");
       ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
      free(value);
    }
    free(key);
    tcbdbcurnext(cursor);
  }
  tcbdbcurdel(cursor);
}

void FragSpectrumScanDatabase::putFSS( ::percolatorInNs::fragSpectrumScan & fss ) {
  assert(bdb);
  *oxdrp << fss;
  xdrrec_endofrecord (&xdr, true);
  ::percolatorInNs::fragSpectrumScan::scanNumber_type key = fss.scanNumber();
  size_t keySize = sizeof(key);
  size_t valueSize(buf.size ());
  if(!tcbdbput(bdb, ( const char * ) &key, keySize, buf.data (), buf.size () ))
  {
    int  errorcode = tcbdbecode(bdb);
    fprintf(stderr, "put error: %s\n", tcbdberrmsg(errorcode));
    exit(EXIT_FAILURE);
  }
  buf.size(0);
  return;
}

extern "C" int
overflow (void* p, char* buf, int n_)
{
  xml_schema::buffer* dst (reinterpret_cast<xml_schema::buffer*> (p));
  size_t n (static_cast<size_t> (n_));
  size_t size (dst->size ());
  size_t capacity (dst->capacity ());
  if (size + n > capacity && size + n < capacity * 2)
    dst->capacity (capacity * 2);
  dst->size (size + n);
  memcpy (dst->data () + size, buf, n);
  return n;
}

extern "C" int
underflow (void* p, char* buf, int n_)
{
  underflow_info* ui (reinterpret_cast<underflow_info*> (p));
  size_t n (static_cast<size_t> (n_));
  size_t size (ui->buf->size () - ui->pos);
  n = size > n ? n : size;
  memcpy (buf, ui->buf->data () + ui->pos, n);
  ui->pos += n;
  return n;
}
