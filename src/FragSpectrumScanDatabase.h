#ifndef FRAGSPECTRUMSCANDATABASE_H
#define FRAGSPECTRUMSCANDATABASE_H

#include <memory>
#include <tcbdb.h>
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "percolator_in.hxx"

class serializer;


class FragSpectrumScanDatabase {
public:
FragSpectrumScanDatabase();
 ~FragSpectrumScanDatabase(){} 

bool init( std::string filename );
std::auto_ptr< ::percolatorInNs::fragSpectrumScan> getFSS( unsigned int scanNr ); 
 std::auto_ptr< ::percolatorInNs::fragSpectrumScan> deserializeFSSfromBinary( char * value, int valueSize );
 void putFSS( ::percolatorInNs::fragSpectrumScan & fss );
 void savePsm( unsigned int scanNr, double observedMassCharge, std::auto_ptr< percolatorInNs::peptideSpectrumMatch > psm_p );
 void print( serializer & ser );
protected:

  XDR xdr;
  xml_schema::buffer buf;

  TCBDB *bdb;

  // is scoped_ptr possible here?
  std::auto_ptr< xml_schema::ostream<XDR> > oxdrp;




};

#endif
