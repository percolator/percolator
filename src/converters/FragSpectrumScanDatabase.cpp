#include "FragSpectrumScanDatabase.h"

#ifndef __BOOSTDB__
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

#endif

FragSpectrumScanDatabase::FragSpectrumScanDatabase(string id_par) :
    bdb(0),scan2rt(0) {
#ifndef __BOOSTDB__
  xdrrec_create_p xdrrec_create_ = reinterpret_cast<xdrrec_create_p> (::xdrrec_create);
  xdrrec_create_ (&xdr, 0, 0, reinterpret_cast<char*> (&buf), 0, &overflow);
  xdr.x_op = XDR_ENCODE;
  std::auto_ptr< xml_schema::ostream<XDR> > tmpPtr(new xml_schema::ostream<XDR>(xdr)) ;
  assert(tmpPtr.get());
  oxdrp=tmpPtr;
  if(id_par.empty()) id = "no_id"; else id = id_par;
#else
  bdb = new mapdb();
#endif
}


void FragSpectrumScanDatabase::savePsm( unsigned int scanNr,
    std::auto_ptr< percolatorInNs::peptideSpectrumMatch > psm_p ) {

  std::auto_ptr< ::percolatorInNs::fragSpectrumScan>  fss = getFSS( scanNr );
  // if FragSpectrumScan does not yet exist, create it
  if ( ! fss.get() ) {
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan>
    fs_p( new ::percolatorInNs::fragSpectrumScan(scanNr));
    fss = fs_p;
  }
  // add the psm to the FragSpectrumScan
  fss->peptideSpectrumMatch().push_back(psm_p);
  putFSS( *fss );
  return;
}

//   The function "tcbdbopen" in Tokyo Cabinet does not have O_EXCL as is
//   possible in the unix system call open (see "man 2 open"). This may be a
//   security issue if the filename to the Tokyo cabinet database is in a
//   directory that other users have write access to. They could add a symbolic
//   link pointing somewhere else. It would be better if Tokyo Cabinet would
//   fail if the database existed in our case when we use a temporary file.
bool FragSpectrumScanDatabase::init(std::string fileName) {

  #if defined __LEVELDB__
    options.create_if_missing = true;
    options.error_if_exists = true;
    options.max_open_files = 100;
    options.write_buffer_size = 4194304*2; //8 MB
    options.block_size = 4096*4; //16K
    
    leveldb::Status status = leveldb::DB::Open(options, fileName.c_str(), &bdb);
    if (!status.ok()){ 
      std::cerr << status.ToString() << endl;
      exit(EXIT_FAILURE);
    }
    bool ret = status.ok();
  #elifdef __TOKYODB__
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
  #else
    bool ret = true;
  #endif
    
  return ret;
}

void FragSpectrumScanDatabase::terminte(){
  #if defined __LEVELDB__
    delete(bdb);
  #elifdef __TOKYODB__
    tcbdbdel(bdb);   
  #endif
}

bool FragSpectrumScanDatabase::initRTime(map<int, vector<double> >* scan2rt_par) {
  // add pointer to retention times table in sqt2pin (if any)
  scan2rt=scan2rt_par;
}

#ifndef __BOOSTDB__
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
  //TOFIX it gives too many arguments in MINGW
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
#endif

std::auto_ptr< ::percolatorInNs::fragSpectrumScan> FragSpectrumScanDatabase::getFSS( unsigned int scanNr ) {
  #if defined __LEVELDB__
    assert(bdb);
    std::string skey = boost::lexical_cast<std::string>(scanNr);
    leveldb::Slice s1(skey);
    leveldb::Iterator* itr = bdb->NewIterator(leveldb::ReadOptions());
    itr->Seek(s1);
    if(!itr->Valid() || s1 != itr->key()){
      delete itr;
      return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
    }
    char *retvalue = const_cast<char*>(itr->value().data());
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret(deserializeFSSfromBinary(retvalue,itr->value().size()));
    delete itr;
    return ret;
  #elifdef __TOKYODB__
    assert(bdb);
    int valueSize = 0;
    char * value = ( char * ) tcbdbget(bdb, ( const char* ) &scanNr, sizeof( scanNr ), &valueSize);
    if(!value) {
      return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
    }
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret(deserializeFSSfromBinary(value,valueSize));  
    free(value);
    return ret;
  #else
    assert(bdb);
    mapdb::iterator it;
    it = bdb->find(scanNr);
    if(it == bdb->end()){
      return std::auto_ptr< ::percolatorInNs::fragSpectrumScan> (NULL);
    }       
    std::istringstream istr (it->second);
    text_iarchive ia (istr);
    xml_schema::istream<text_iarchive> is (ia);
    std::auto_ptr< ::percolatorInNs::fragSpectrumScan> ret (new ::percolatorInNs::fragSpectrumScan (is));    
    return ret;
  #endif
  
}

void FragSpectrumScanDatabase::print(serializer & ser) {
  #if defined __LEVELDB__
    assert(bdb);
    leveldb::Iterator* it = bdb->NewIterator(leveldb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
      char *retvalue = const_cast<char*>(it->value().data());
      std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss(deserializeFSSfromBinary(retvalue,it->value().size()));
      ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
    }
    delete it;
  #elifdef __TOKYODB__
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
	ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
	free(value);
      }
      free(key);
      tcbdbcurnext(cursor);
    }
    tcbdbcurdel(cursor);
  #else
    assert(bdb);
    mapdb::iterator it;
    for (it = bdb->begin(); it != bdb->end(); it++) {
       std::istringstream istr (it->second);
       text_iarchive ia (istr);
       xml_schema::istream<text_iarchive> is (ia);
       std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss (new ::percolatorInNs::fragSpectrumScan (is));  
       ser.next ( PERCOLATOR_IN_NAMESPACE, "fragSpectrumScan", *fss);
    }
  #endif
}

void FragSpectrumScanDatabase::putFSS( ::percolatorInNs::fragSpectrumScan & fss ) {
  #if defined __LEVELDB__
    assert(bdb);
    *oxdrp << fss;
    xdrrec_endofrecord (&xdr, true);
    leveldb::WriteOptions write_options;
    ::percolatorInNs::fragSpectrumScan::scanNumber_type key = fss.scanNumber();
    std::string skey = boost::lexical_cast<std::string>(key);
    leveldb::Slice s2(buf.data(),buf.size());
    leveldb::Slice s1(skey);
    leveldb::Status status = bdb->Put(write_options,s1,s2);
    if(!status.ok()){
      std::cerr << status.ToString() << endl;
      exit(EXIT_FAILURE);
    }
    buf.size(0);
  #elifdef __TOKYODB__
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
  #else
    std::ostringstream ostr;
    text_oarchive oa (ostr);
    xml_schema::ostream<text_oarchive> os (oa);
    os << fss;
    bdb->insert(make_pair<unsigned int, std::string >((unsigned int)fss.scanNumber(),ostr.str()));
  #endif
}

#ifndef __BOOSTDB__
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
#endif