#include "FragSpectrumScanDatabase.h"
#include <MSToolkitTypes.h>
 

FragSpectrumScanDatabase::FragSpectrumScanDatabase(string id_par) :
    scan2rt(0) 
{
  if(id_par.empty()) id = "no_id"; else id = id_par;
}

void FragSpectrumScanDatabase::savePsm( unsigned int scanNr,
    std::auto_ptr< percolatorInNs::peptideSpectrumMatch > psm_p ) 
{
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

bool FragSpectrumScanDatabase::initRTime(map<int, vector<double> >* scan2rt_par) {
  // add pointer to retention times table (if any)
  scan2rt=scan2rt_par;
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
