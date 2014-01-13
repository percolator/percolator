#include "FragSpectrumScanDatabase.h"
//#include <MSToolkitTypes.h>
 

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

void FragSpectrumScanDatabase::printTabFss(std::auto_ptr< ::percolatorInNs::fragSpectrumScan> fss, ofstream &tabOutputStream) {
  int label = 0;
  for (const auto &psm : fss->peptideSpectrumMatch()) {
    if (psm.isDecoy()) {
      label = -1;
    } else {
      label = 1;
    }

    tabOutputStream << psm.id() << '\t' << label << '\t' << fss->scanNumber();
    for (const auto & feature : psm.features().feature()) {
      tabOutputStream << '\t' << feature;
    }
    bool isFirst = true;
    for( const auto & oc : psm.occurence() ) {
      // adding n-term and c-term residues to peptide
      //NOTE the residues for the peptide in the PSMs are always the same for every protein
      if (isFirst) {
        tabOutputStream << '\t' << oc.flankN() << "." << decoratePeptide(psm.peptide()) << "." << oc.flankC();
        isFirst = false;
      }
      tabOutputStream << '\t' << oc.proteinId();
    }
    tabOutputStream << std::endl;
  }
}

std::string FragSpectrumScanDatabase::decoratePeptide(const ::percolatorInNs::peptideType& peptide) {
  std::list<std::pair<int,std::string> > mods;
  std::string peptideSeq = peptide.peptideSequence();
  for(const auto &mod_ref : peptide.modification()){
    std::stringstream ss;
    if (mod_ref.uniMod().present()) {
      ss << "[UNIMOD:" << mod_ref.uniMod().get().accession() << "]";
      mods.push_back(std::pair<int,std::string>(mod_ref.location(),ss.str()));
    }
    if (mod_ref.freeMod().present()) {
      ss << "[" << mod_ref.freeMod().get().moniker() << "]";
      mods.push_back(std::pair<int,std::string>(mod_ref.location(),ss.str()));
    }
  }
  mods.sort(greater<std::pair<int,std::string> >());
  std::list<std::pair<int,std::string> >::const_iterator it;
  for(it=mods.begin();it!=mods.end();++it) {
    peptideSeq.insert(it->first,it->second);
  }
  return peptideSeq;
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
