#ifndef _MSREADER_H
#define _MSREADER_H

#include "Spectrum.h"
#include "MSObject.h"
#include <cstring>
#include "zlib.h"
#include "ramp.h"
#include "ramp_base64.h"
#include <algorithm>

#ifndef _NOSQLITE
#include <sqlite3.h>
#endif


//Macros for 64-bit file support
#ifdef _MSC_VER
//extern "C" int __cdecl _fseeki64(FILE *, __int64, int);
//extern "C" __int64 __cdecl _ftelli64(FILE *);
typedef __int64 f_off;
#define fseek(h,p,o) _fseeki64(h,p,o)
#define ftell(h) _ftelli64(h)
#else

#ifndef _LARGEFILE_SOURCE
#error "need to define _LARGEFILE_SOURCE!!"
#endif    /* end _LARGEFILE_SOURCE */

typedef off_t f_off;
#define fseek(h,p,o) fseeko(h,p,o)
#define ftell(h) ftello(h)

#endif /* end _MSC_VER */ 


using namespace std;

class MSReader {
 public:
  //Constructors & Destructors
  MSReader();
  ~MSReader();

  //Functions
  void appendFile(char* c, bool text, Spectrum& s);
  void appendFile(char* c, bool text, MSObject& m);
  void appendFile(char* c, Spectrum& s);
  void appendFile(char* c, MSObject& m);
  
  MSFileFormat checkFileFormat(char *fn);

  MSHeader& getHeader();
  //Spectrum readBinaryFile(char* c, Spectrum& s, int scNum=0);
  //Spectrum readMSFile(char* c,int scNum=0);
	
  MSSpectrumType getFileType();
  int getPercent();
  int getLastScan();

  //get specific header informations
  void setPrecision(int i, int j);
  void setPrecisionInt(int i);
  void setPrecisionMZ(int i);
  void writeFile(char* c, bool text, MSObject& m);
  void writeFile(char* c, MSFileFormat ff, MSObject& m, char* sha1Report='\0');

  bool readFile(char* c, bool text, Spectrum& s, int scNum=0);
  bool readFile(char* c, MSFileFormat f, Spectrum& s, int scNum=0);
  bool readFile(char* c, Spectrum& s, int scNum=0);
  
  void setFilter(vector<MSSpectrumType>& m);
  void setFilter(const MSSpectrumType& m);

  //File compression
  void setCompression(bool b);

  //for Sqlite
  void createIndex(); 

 protected:

 private:
  //Data Members
  FILE *fileIn;
  MSHeader header;
  int headerIndex;
  MSSpectrumType fileType;
  f_off lEnd;
  f_off lPivot;
  f_off lFWidth;
  int iIntensityPrecision;
  int iMZPrecision;
  int iVersion;
  int iFType;
  MSFileFormat lastFileFormat;

  //File compression
  bool compressMe;

  //mzXML support variables;
  ramp_fileoffset_t  *pScanIndex;
  RAMPFILE  *rampFileIn;
  bool rampFileOpen;
  int rampLastScan;
  int rampIndex;
  vector<MSSpectrumType> filter;

  //Functions
  void closeFile();
  int openFile(char* c, bool text=false);
  bool findSpectrum(int i);
  void readCompressSpec(FILE* fileIn, MSScanInfo& ms, Spectrum& s);
  void readSpecHeader(FILE* fileIn, MSScanInfo& ms);

  void writeBinarySpec(FILE* fileOut, Spectrum& s);
  void writeCompressSpec(FILE* fileOut, Spectrum& s);
  void writeTextSpec(FILE* fileOut, Spectrum& s);
  void writeSpecHeader(FILE* fileOut, bool text, Spectrum& s);
  
  //support for sqlite
  #ifndef _NOSQLITE
  bool readSqlite(char* c, Spectrum& s, int scNum);
  void getUncompressedPeaks(Spectrum& s, int& numPeaks, int& mzLen, Byte* comprM, int& intensityLen, Byte* comprI);
  int curIndex;  //remember where we are
  int lastScanNumber;
  int lastIndex;
  sqlite3* db;
  void sql_stmt(const char* stmt);
  bool executeSqlStmt(Spectrum& s, char* zSql);
  void appendFile(Spectrum& s);
  void writeSqlite(char* c, MSObject& m, char* sha1Report);
  void readChargeTable(int scanID, Spectrum& s);
  vector<int> estimateCharge(Spectrum& s);
  #endif 

};

#endif

