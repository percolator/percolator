#include "MSReader.h"
#include <iostream>
using namespace std;

MSReader::MSReader(){
  fileIn=NULL;
  rampFileIn=NULL;
  iIntensityPrecision=1;
  iMZPrecision=4;
  //filter=Unspecified;
  rampFileOpen=false;
  compressMe=false;
  rawFileOpen=false;
  exportMGF=false;
  highResMGF=false;
  iFType=0;
  iVersion=0;
  for(int i=0;i<16;i++)	strcpy(header.header[i],"\0");
  headerIndex=0;

  #ifdef _MSC_VER
  CoInitialize( NULL );
	HRESULT hr = m_Raw.CreateInstance("XRawfile.XRawfile.1");
  if (!FAILED(hr)) bRaw=true;
  else bRaw=false;
  rawCurSpec=0;
  rawTotSpec=0;
  rawAvg=false;
  rawAvgWidth=1;
  rawAvgCutoff=1000;
  rawLabel=false;
  rawUserFilterExact=true;
  strcpy(rawUserFilter,"");
  #endif

  #ifndef _NOSQLITE
  db = NULL;
  lastIndex=-1;
  lastScanNumber=-1;
  #endif
}

MSReader::~MSReader(){
  closeFile();
  if(rampFileOpen) {
    rampCloseFile(rampFileIn);
    free(pScanIndex);
  }

  #ifdef _MSC_VER
  if(bRaw){
    if(rawFileOpen) m_Raw->Close();
    m_Raw.Release();
    m_Raw=NULL;
  }
  #endif

  #ifndef _NOSQLITE
  if(db != NULL)  sqlite3_close(db);
  #endif
}

void MSReader::closeFile(){
  if(fileIn!=NULL) fclose(fileIn);
}

MSHeader& MSReader::getHeader(){
  return header;
}

/* 0 = File opened correctly
   1 = Could not open file
*/
int MSReader::openFile(const char *c,bool text){
	int i;

	if(text) fileIn=fopen(c,"rt");
	else fileIn=fopen(c,"rb");

  if(fileIn==NULL) {
		for(i=0;i<16;i++) strcpy(header.header[i],"\0");
    headerIndex=0;
    fileType=Unspecified;
    return 1;
  } else {
    fileType=Unspecified;

		//if we don't have the eof position, get it here.
		fseek(fileIn,0,2);
		lEnd=ftell(fileIn);

		lPivot = 0;
    lFWidth = lEnd/2;

		fseek(fileIn,0,0);

		if(text){
			for(i=0;i<16;i++) strcpy(header.header[i],"\0");
			headerIndex=0;
		} else {
      fread(&iFType,4,1,fileIn);
      fread(&iVersion,4,1,fileIn);
			fread(&header,sizeof(MSHeader),1,fileIn);
		};

	  return 0;
  };
};

MSSpectrumType MSReader::getFileType(){
  return fileType;
}


bool MSReader::readFile(const char *c, bool text, Spectrum& s, int scNum){
  MSScanInfo ms;
  Peak_T p;
  ZState z;
  EZState ez;
  int i;

  //variables for text reading only
  bool firstScan = false;
  bool bScan = true;
  bool bDoneHeader = false;
  char tstr[256];
  char ch;
  char *tok;
  //off_t fpoint;

  //variables for compressed files
  uLong mzLen, intensityLen;

  //clear any spectrum data
  s.clear();

  //check for valid file and if we can access it
  if(c!=NULL){
    closeFile();
    if(openFile(c,text)==1) return false;
    lastFileFormat = checkFileFormat(c);
  } else if(fileIn==NULL) {
    return false;
  }

  //set the filetype
  switch(lastFileFormat){
  case ms2:
  case cms2:
  case bms2:
    s.setFileType(MS2);
    break;
  case zs:
    s.setFileType(ZS);
    break;
  case uzs:
    s.setFileType(UZS);
    break;
  case ms1:
  case cms1:
  case bms1:
    s.setFileType(MS1);
    break;
  default:
    s.setFileType(Unspecified);
    break;
  }

	//Handle binary and text files differently
  if(!text){

    //if binary file, read scan info sequentially, skipping to next scan if requested
    //fread(&ms,sizeof(MSScanInfo),1,fileIn);
    readSpecHeader(fileIn,ms);

    if(scNum!=0){

      fseek(fileIn,sizeof(MSHeader)+8,0);

      //fread(&ms,sizeof(MSScanInfo),1,fileIn);
      readSpecHeader(fileIn,ms);

      while(ms.scanNumber[0]!=scNum){

        fseek(fileIn,ms.numZStates*12,1);
        fseek(fileIn,ms.numEZStates*20,1);

	      if(compressMe){
	        fread(&i,4,1,fileIn);
	        mzLen = (uLong)i;
	        fread(&i,4,1,fileIn);
	        intensityLen = (uLong)i;
	        fseek(fileIn,mzLen+intensityLen,1);
	      } else {
	        fseek(fileIn,ms.numDataPoints*12,1);
	      }

	      //fread(&ms,sizeof(MSScanInfo),1,fileIn);
        readSpecHeader(fileIn,ms);
	      if(feof(fileIn)) return false;
      }
    }
    if(feof(fileIn)) return false;

		//read any charge states (for MS2 files)
    for(i=0;i<ms.numZStates;i++){
      fread(&z.z,4,1,fileIn);
      fread(&z.mz,8,1,fileIn);
      s.addZState(z);
    }

    for(i=0;i<ms.numEZStates;i++){
      fread(&ez.z,4,1,fileIn);
      fread(&ez.mh,8,1,fileIn);
      fread(&ez.pRTime,4,1,fileIn);
      fread(&ez.pArea,4,1,fileIn);
      s.addEZState(ez);
    }

    s.setScanNumber(ms.scanNumber[0]);
    s.setScanNumber(ms.scanNumber[1],true);
    s.setRTime(ms.rTime);
    s.setMZ(ms.mz);
    s.setBPI(ms.BPI);
    s.setBPM(ms.BPM);
    s.setConversionA(ms.convA);
    s.setConversionB(ms.convB);
    s.setIonInjectionTime(ms.IIT);
    s.setTIC(ms.TIC);

    //read compressed data to the spectrum object
    if(compressMe) {

      readCompressSpec(fileIn,ms,s);

      //or read binary data to the spectrum object
    } else {
      for(i=0;i<ms.numDataPoints;i++){
	      fread(&p.mz,8,1,fileIn);
	      fread(&p.intensity,4,1,fileIn);
	      //cout << p.mz << " " << p.intensity << endl;
	      s.add(p);
      }
    }

    //return success
    return true;

  } else {

    //if reading text files, some parsing is required.
    while(true){

      //stop when you reach the end of the file
      if(feof(fileIn)) {

        //Special case: when doing binary search, end of file might mean to search
        //the othere end of the file.
        if(scNum != 0){
	        if(s.getScanNumber() != scNum) {
	          bScan=findSpectrum(-1);
            s.clear();
	          s.setScanNumber(0);
	          if(bScan==false) return false;
          } else {
            break;
          }
	      } else {
          break;
        }
      }

      //scan next character in the file
      ch=fgetc(fileIn);
      ungetc(ch,fileIn);

      switch(ch){
      case 'D':
	      //D lines are ignored
	      fgets(tstr,256,fileIn);
	      break;

      case 'H':
	      //Header lines are recorded as strings up to 16 lines at 256 characters each
	      fgets(tstr,256,fileIn);
	      if(!bDoneHeader) {
	        tok=strtok(tstr," \t\n\r");
	        tok=strtok(NULL,"\n\r");
	        strcat(tok,"\n");
	        if(headerIndex<16) strcpy(header.header[headerIndex++],tok);
	        else cout << "Header too big!!" << endl;
	      }
	      break;

      case 'I':
	      //I lines are recorded only if they contain retention times
        fgets(tstr,256,fileIn);
        tok=strtok(tstr," \t\n\r");
        tok=strtok(NULL," \t\n\r");
        if(strcmp(tok,"RTime")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setRTime((float)atof(tok));
        }	else if(strcmp(tok,"TIC")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setTIC((float)atof(tok));
        }	else if(strcmp(tok,"IIT")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setIonInjectionTime((float)atof(tok));
        }	else if(strcmp(tok,"BPI")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setBPI((float)atof(tok));
        }	else if(strcmp(tok,"BPM")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setBPM((float)atof(tok));
        }	else if(strcmp(tok,"ConvA")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setConversionA(atof(tok));
        }	else if(strcmp(tok,"ConvB")==0) {
          tok=strtok(NULL," \t\n\r,");
          s.setConversionB(atof(tok));
        } else if(strcmp(tok,"EZ")==0) {
          tok=strtok(NULL," \t\n\r,");
          ez.z=atoi(tok);
          tok=strtok(NULL," \t\n\r,");
          ez.mh=atof(tok);
          tok=strtok(NULL," \t\n\r,");
          ez.pRTime=(float)atof(tok);
          tok=strtok(NULL," \t\n\r,");
          ez.pArea=(float)atof(tok);
          s.addEZState(ez);
        }
        break;

      case 'S':
	      //Scan numbers are recorded and mark all following data is spectrum data
	      //until the next tag

	      //Reaching an S tag also indicates there are no more header lines
	      bDoneHeader=true;

	      if(firstScan) {
	        //if we are here, a desired scan was read and we just reached the next scan tag
	        //therefore, stop reading further.
	        return true;

	      } else {
	        fgets(tstr,256,fileIn);
	        tok=strtok(tstr," \t\n\r");
	        tok=strtok(NULL," \t\n\r");
          s.setScanNumber(atoi(tok));
	        tok=strtok(NULL," \t\n\r");
	        s.setScanNumber(atoi(tok),true);
	        tok=strtok(NULL," \t\n\r");
	        if(tok!=NULL)	s.setMZ(atof(tok));
	        if(scNum != 0){
	          if(s.getScanNumber() != scNum) {
	            if(s.getScanNumber()<scNum) bScan=findSpectrum(1);
	            else bScan=findSpectrum(-1);
              s.clear();
	            s.setScanNumber(0);
	            if(bScan==false) return false;
	            break;
	          }
	        }
	        firstScan=true;
	      }
	      break;

      case 'Z':
	      //Z lines are recorded for MS2 files

        //don't record z-lines unless this is a scan number that is wanted
        if(!firstScan){
	  fgets(tstr,256,fileIn);
 	  break;
	}
	      fgets(tstr,256,fileIn);
	      tok=strtok(tstr," \t\n\r");
	      tok=strtok(NULL," \t\n\r");
	      z.z=atoi(tok);
	      tok=strtok(NULL," \t\n\r");
	      z.mz=atof(tok);
	      s.addZState(z);
	      break;

      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
	      //lines beginning with numbers are data; if they belong to a scan we are not
	      //interested in, we ignore them.
	      if(scNum != 0){
	        if(s.getScanNumber()!=scNum) {
	          fgets(tstr,256,fileIn);
	          break;
	        }
	      }
	      //otherwise, read in the line
	      fscanf(fileIn,"%lf %f\n",&p.mz,&p.intensity);
	      s.add(p);
	      break;

      default:
	      //if the character is not recognized, ignore the entire line.
        fgets(tstr,256,fileIn);
	      //fscanf(fileIn,"%s\n",tstr);
	      break;
      }
    }

  }

  return true;

}


bool MSReader::findSpectrum(int i){

  if(i==0){
    lPivot = lEnd/2;
    lFWidth = lPivot/2;
  } else if(i==-1){
    lPivot -= lFWidth;
    lFWidth /= 2;
  } else {
    lPivot += lFWidth;
    lFWidth /= 2;
  }

  fseek(fileIn,lPivot,0);
  return (lFWidth>0 && lPivot>0 && lPivot<lEnd);

}

int MSReader::getLastScan(){
  if(rampFileIn!=NULL){
    return (rampLastScan);
  }
  #ifndef _NOSQLITE
  if(db != 0){
    return lastScanNumber;
  }
  #endif
  return -1;
}

int MSReader::getPercent(){
  if(fileIn!=NULL){
    return (int)((double)ftell(fileIn)/lEnd*100);
  }
  if(rampFileIn!=NULL){
    return (int)((double)rampIndex/rampLastScan*100);
  }
#ifdef _MSC_VER
  if(rawFileOpen){
    return (int)((double)rawCurSpec/rawTotSpec*100);
  }
#endif
  return -1;
}

void MSReader::writeFile(const char* c, bool text, MSObject& m){

  FILE* fileOut;
  int i;

  //if a filename isn't specified, check to see if the
  //MSObject has a filename.
  if(c == NULL) {
    return;
  } else {
    if(text) fileOut=fopen(c,"wt");
    else fileOut=fopen(c,"wb");
  }

  //output file header lines;
  if(text){
    if(exportMGF){
      //MGF file header is here
      fprintf(fileOut,"COM=Generated in the MSToolkit\n");
      if(!highResMGF) fprintf(fileOut,"CHARGE=2+ and 3+\n");
    } else {
      //MSx/SQT file header is here
      for(i=0;i<16;i++){
        if(m.getHeader().header[i][0]!='\0') {
          fputs("H\t",fileOut);
          fputs(m.getHeader().header[i],fileOut);
        }
      }
    }
  } else {
    //version 0 or 1 has basic stats
    //version 2 adds BPI,BPM,TIC,IIT,ConvA,ConvB
    //version 3 adds EZ lines
    fwrite(&iFType,4,1,fileOut); //file type
    i=3;
    fwrite(&i,4,1,fileOut); //version number - in case we change formats
		fwrite(&m.getHeader(),sizeof(MSHeader),1,fileOut);
	}

	//output spectra;
  for(i=0;i<m.size();i++){

		//output spectrum header
		writeSpecHeader(fileOut,text,m.at(i));

		//output scan
		if(text){
			writeTextSpec(fileOut,m.at(i));
		} else if(compressMe){
			writeCompressSpec(fileOut,m.at(i));
		} else {
			writeBinarySpec(fileOut,m.at(i));
		}

  }

	fclose(fileOut);
}

void MSReader::writeFile(const char* c, MSFileFormat ff, MSObject& m, char* sha1Report){

  switch(ff){
  case mgf:
    exportMGF=true;
    setCompression(false);
    writeFile(c,true,m);
    exportMGF=false;
    break;
  case ms1:
  case ms2:
  case  zs:
  case uzs:
    exportMGF=false;
    setCompression(false);
    writeFile(c,true,m);
    break;
  case psm:
    #ifndef _NOSQLITE
    writeSqlite(c,m, sha1Report);
    #endif
    break;
  case mzXML:
  case mzData:
    cout << "Cannot write mzXML or mzData formats. Nothing written." << endl;
    break;
  case bms1:
    exportMGF=false;
    setCompression(false);
    iFType=1;
    writeFile(c,false,m);
    break;
  case bms2:
    exportMGF=false;
    setCompression(false);
    iFType=3;
    writeFile(c,false,m);
    break;
  case cms1:
    exportMGF=false;
    setCompression(true);
    iFType=2;
    writeFile(c,false,m);
    break;
  case cms2:
    exportMGF=false;
    setCompression(true);
    iFType=4;
    writeFile(c,false,m);
    break;
  default:
    cout << "Unknown file format. Nothing written." << endl;
    break;
  }

}

#ifndef _NOSQLITE
void MSReader::writeSqlite(const char* c, MSObject& m, char* sha1Report)
{

  //open the database for write
  sqlite3_open(c, &db);
  if(db == 0)
    {
      cout<<"Error open database "<<c<<endl;
      return;
    }

  sql_stmt("PRAGMA synchronous=OFF");
  sql_stmt("PRAGMA cache_size=750000");
  sql_stmt("PRAGMA temp_store=MEMORY");

  //create two tables (msRun and msScan)
  char zSql[8192];

  strcpy(zSql, "create table msRun(id INTEGER primary key autoincrement not null,"
	 "filename VARCHAR(150), sha1Sum VARCHAR(100), creationTime VARCHAR(255), extractor VARCHAR(255),"
	 "extractorVersion VARCHAR(100), instrumentType VARCHAR(100), instrumentVendor TEXT, instrumentSN TEXT,"
	 "acquisitionMethod TEXT, originalFileType TEXT, separateDigestion CHAR(1), uploadDate TEXT, comment TEXT)");
  sql_stmt(zSql);
  zSql[0]='\0';

  strcpy(zSql,"create table msScan(id INTEGER primary key autoincrement not null,"
	 "runID INTEGER, startScanNumber INTEGER, endScanNumber INTEGER, level INTEGER, precursorMZ REAL, precursorCharge INTEGER,"
	 "preScanID INTEGER, preScanNumber INTEGER, retentionTime REAL, fragmentationType VARCHAR(100),isCentroid CHAR(1), peakCount INTEGER)");
  sql_stmt(zSql);
  zSql[0]='\0';

  strcpy(zSql, "create table msScanData(scanID INTEGER, peakMZ BLOB, peakIntensity BLOB)");
  sql_stmt(zSql);
  zSql[0]='\0';

  strcpy(zSql,"create table MS2FileScanCharge(id INTEGER primary key autoincrement not null,"
	 "scanID INTEGER, charge INTEGER, mass REAL)");
  sql_stmt(zSql);
  zSql[0]='\0';

  //get the creationTime, sha1Sum etc.
  string fileCreateTime="=";
  string instrumentType="=";
  for(int i=0; i<16; i++)
    {
      if(m.getHeader().header[i] != '\0')
	{
	  string headerLine = m.getHeader().header[i];
	  if(headerLine.find("CreationDate") != string::npos)
	    {
	      size_t pos = headerLine.find('\t');
	      fileCreateTime = headerLine.substr(pos+1);
	    }
	  if(headerLine.find("InstrumentType") != string::npos)
	    {
	      size_t pos = headerLine.find('\t');
	      instrumentType = headerLine.substr(pos+1);
	    }
	}
    }


  //insert into table msRun
  sprintf(zSql,"insert into msRun(filename, sha1Sum, creationTime, extractor,extractorVersion, instrumentType) values('%s','%s','%s','%s','%s','%s')",
	  c,
	  sha1Report,
	  fileCreateTime.c_str(),
	  "MakeMS2",
	  "1.0",
	  instrumentType.c_str());

  sql_stmt(zSql);

}
#endif

void MSReader::appendFile(char* c, bool text, Spectrum& s){
  FILE* fileOut;

  if(c == NULL) return;

  if(text)fileOut=fopen(c,"at");
  else fileOut=fopen(c,"ab");

  //output spectrum header
  writeSpecHeader(fileOut,text,s);

  //output spectrum
  if(text){
    writeTextSpec(fileOut,s);
  } else if(compressMe){
    writeCompressSpec(fileOut,s);
  } else {
    writeBinarySpec(fileOut,s);
  }

  fclose(fileOut);

}

void MSReader::appendFile(char* c, Spectrum& s){
  MSFileFormat ff;
  FILE* fileOut;

  if(c == NULL) return;
  ff=checkFileFormat(c);

  switch(ff){
  case mgf:
    exportMGF=true;
    fileOut=fopen(c,"at");
    writeTextSpec(fileOut,s);
    fclose(fileOut);
    exportMGF=false;
    break;
  case ms1:
  case ms2:
  case  zs:
  case uzs:
    fileOut=fopen(c,"at");
    writeSpecHeader(fileOut,true,s);
    writeTextSpec(fileOut,s);
    fclose(fileOut);
    break;
  case bms1:
  case bms2:
    fileOut=fopen(c,"ab");
    writeSpecHeader(fileOut,false,s);
    writeBinarySpec(fileOut,s);
    fclose(fileOut);
    break;
  case cms1:
  case cms2:
    fileOut=fopen(c,"ab");
    writeSpecHeader(fileOut,false,s);
    writeCompressSpec(fileOut,s);
    fclose(fileOut);
    break;
  case psm:
    #ifndef _NOSQLITE
    appendFile(s);
    #endif
    break;
  default:
    cout << "Cannot append file: unknown or unsupported file type." << endl;
    break;
  }

}

//private function for insert a scan into msScan table
#ifndef _NOSQLITE
void MSReader::appendFile(Spectrum& s)
{

  if(db == 0)
    return;

  int j;

  //file compression
  int err;
  uLong len;
  unsigned char *comprM, *comprI;
  uLong comprLenM, comprLenI;
  double *pD;
  float *pF;
  uLong sizeM;
  uLong sizeI;

  //Build arrays to hold scan prior to compression

  pD = new double[s.size()];
  pF = new float[s.size()];
  for(j=0;j<s.size();j++){
    pD[j]=s.at(j).mz;
    pF[j]=s.at(j).intensity;
  }

  //compress mz
  len = (uLong)s.size()*sizeof(double);
  sizeM = len;
  comprLenM = compressBound(len);
  comprM = (unsigned char*)calloc((uInt)comprLenM, 1);
  err = compress(comprM, &comprLenM, (const Bytef*)pD, len);

  //compress intensity
  len = (uLong)s.size()*sizeof(float);
  sizeI = len;
  comprLenI = compressBound(len);
  comprI = (unsigned char*)calloc((uInt)comprLenI, 1);
  err = compress(comprI, &comprLenI, (const Bytef*)pF, len);

   //insert into table
  char zSql[8192];
  int charge;

  vector<int> chgs;
  if(s.getMsLevel() >1)
    {
      chgs = estimateCharge(s);
    }

  if(s.sizeZ() == 1)
    charge = s.atZ(0).z;
  else if(s.sizeZ() > 1)
    charge = 0;
  else
    {
      if(chgs.size() == 1)
	charge = chgs.at(0);
      else
	charge = 0;
    }


  MSActivation act = s.getActivationMethod();
  string actMethod;
  switch(act){
  case ETD:
    actMethod="ETD";
    break;
  case CID:
    actMethod="CID";
    break;
  case ECD:
    actMethod="ECD";
    break;
  case PQD:
    actMethod = "PQD";
    break;
  case HCD:
    actMethod = "HCD";
    break;
  case na:
    actMethod="UNKNOWN";
    break;
  default:
    actMethod="UNKNOWN";
    break;
  }


  sprintf(zSql, "insert into msScan(runID,startScanNumber,endScanNumber,level,precursorMZ, precursorCharge,retentionTime,fragmentationType,peakCount) "
	  "values (1,%d, %d,%d, %f, %d, %f,'%s', %d)",
          s.getScanNumber(),
	  s.getScanNumber(true),
	  s.getMsLevel(),
          s.getMZ(),
          charge,
          s.getRTime(),
	  actMethod.c_str(),
          s.size());


  //cout<<zSql<<endl;
  sql_stmt(zSql);
  zSql[0]='\0';

  //get scanID
  strcpy(zSql, "select MAX(id) from msScan");
  int rc,iRow, iCol;
  char** result;

  int lastScanID;
   rc = sqlite3_get_table(db, zSql, &result, &iRow, &iCol, 0);

   if(rc == SQLITE_OK)
     {
       lastScanID=atoi(result[1]);

     }
   else
     {
       cout<<"Can't execute the SQL statement"<<zSql<<endl;
     }

   zSql[0]='\0';

   //insert into msScanData
   sprintf(zSql, "insert into msScanData values(%d, ?, ?)",
	   lastScanID);

  sqlite3_stmt *pStmt;


  rc = sqlite3_prepare(db, zSql, -1, &pStmt, 0);
  if( rc!=SQLITE_OK ){
    cout<<"can't prepare SQL statement!"<<rc<<endl;
    exit(1);
  }


  sqlite3_bind_blob(pStmt, 1, comprM, (int)comprLenM, SQLITE_STATIC);
  sqlite3_bind_blob(pStmt, 2, comprI, (int)comprLenI, SQLITE_STATIC);

  rc = sqlite3_step(pStmt);
  rc = sqlite3_finalize(pStmt);

  free(comprM);
  free(comprI);
  delete [] pD;
  delete [] pF;

  zSql[0]='\0';
  //insert into MS2FileScanCharge
  int chg;
  double MH;
  if(s.getMsLevel() > 1)
    {
      if(s.sizeZ() > 0)
	{
	  for(int i=0; i<s.sizeZ(); i++)
	    {
	      chg = s.atZ(i).z;
	      MH = s.getMZ()*chg-(chg-1)*1.008;
	      sprintf(zSql, "insert into MS2FileScanCharge(scanID, charge, mass) values(%d, %d, %f)",
		     lastScanID,
		     chg,
		     MH);
	      sql_stmt(zSql);
	    }

	}
      else
	{

	  for(int i=0; i<chgs.size(); i++)
	    {
	      chg=chgs.at(i);
	      MH = s.getMZ()*chg -(chg-1)*1.008;
	      sprintf(zSql, "insert into MS2FileScanCharge(scanID, charge, mass) values (%d, %d, %f)",
		     lastScanID,
		     chg,
		     MH);
	      sql_stmt(zSql);
	    }
	}
    }
}

vector<int> MSReader::estimateCharge(Spectrum& s)
{
  vector<int> chgs;
  float totalIntensity = s.getTotalIntensity();
  float beforeMZIntensity=0;

  double preMZ = s.getMZ();

  for(int i=0; i<s.size(); i++)
    {
      if(s.at(i).mz <= preMZ)
	beforeMZIntensity += s.at(i).intensity;
      else
	break;
    }

  if(beforeMZIntensity/totalIntensity >= 0.95)
    {
      chgs.push_back(1);
    }
  else
    {
      chgs.push_back(2);
      chgs.push_back(3);
    }
  return chgs;
}


void MSReader::createIndex()
{
  //create index for msScan table
  char* stmt1 = "create index idxScanNumber on msScan(startScanNumber)";
  sql_stmt(stmt1);

}
#endif

void MSReader::appendFile(char* c, bool text, MSObject& m){

  FILE* fileOut;
  int i;

  //if a filename isn't specified, check to see if the
  //MSObject has a filename.
  if(c == NULL) {
		return;
  } else {
		if(text) fileOut=fopen(c,"at");
		else fileOut=fopen(c,"ab");
	};

  //output spectra;
  for(i=0;i<m.size();i++){

		//output spectrum header
		writeSpecHeader(fileOut,text,m.at(i));

		//output spectrum
		if(text){
			writeTextSpec(fileOut,m.at(i));
		} else if(compressMe){
			writeCompressSpec(fileOut,m.at(i));
		} else {
			writeBinarySpec(fileOut,m.at(i));
		};

	};

	fclose(fileOut);
};

void MSReader::appendFile(char* c, MSObject& m){

  MSFileFormat ff;
  FILE* fileOut;
  int i;

  if(c == NULL) return;
  ff=checkFileFormat(c);

  switch(ff){
    case mgf:
      exportMGF=true;
      fileOut=fopen(c,"at");
      for(i=0;i<m.size();i++) writeTextSpec(fileOut,m.at(i));
      fclose(fileOut);
      exportMGF=false;
      break;
    case ms1:
    case ms2:
    case  zs:
    case uzs:
	    fileOut=fopen(c,"at");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,true,m.at(i));
        writeTextSpec(fileOut,m.at(i));
      }
      fclose(fileOut);
      break;
    case bms1:
    case bms2:
      fileOut=fopen(c,"ab");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,false,m.at(i));
        writeBinarySpec(fileOut,m.at(i));
      }
      fclose(fileOut);
      break;
    case cms1:
    case cms2:
      fileOut=fopen(c,"ab");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,false,m.at(i));
        writeCompressSpec(fileOut,m.at(i));
      }
      fclose(fileOut);
      break;
    default:
      cout << "Cannot append file: unknown or unsupported file type." << endl;
      break;
  }

}

void MSReader::setPrecision(int i, int j){
  iIntensityPrecision=i;
  iMZPrecision=j;
}

void MSReader::setPrecisionInt(int i){
  iIntensityPrecision=i;
}

void MSReader::setPrecisionMZ(int i){
  iMZPrecision=i;
}

bool MSReader::readFile(const char* c, Spectrum& s, int scNum){

  if(c!=NULL) lastFileFormat = checkFileFormat(c);
  return readFile(c,lastFileFormat,s,scNum);

}

#ifndef _NOSQLITE
bool MSReader::readSqlite(const char* c, Spectrum& s, int scNum)
{

  if(c != NULL)
    {
      sqlite3_open(c, &db);
      if(db == 0)
	{
	  cout<<"Error open database "<<c<<endl;
	  return false;
	}

      sql_stmt("PRAGMA synchronous=OFF");
      sql_stmt("PRAGMA cache_size=750000");
      sql_stmt("PRAGMA temp_store=MEMORY");

      char zSql[1024];
      strcpy(zSql, "select MAX(id),MAX(startScanNumber) from msScan");
      int iRow, iCol, rc;
      char** result;


      rc = sqlite3_get_table(db, zSql, &result, &iRow, &iCol, 0);

      if(rc == SQLITE_OK)
        {
	  lastIndex = atoi(result[2]);
          lastScanNumber=atoi(result[3]);
	}

      curIndex = 0;
    }

  s.clear();

  char zSql[2048];
  int rc;


  if(scNum != 0)
    {
      //first get the curIndex

      if(scNum > lastScanNumber)
	{
	  cout<<"Specified scan number doesn't exist!"<<endl;
	  return false;
	}

      sprintf(zSql, "select id from msScan where startScanNumber=%d", scNum);
      int iRow, iCol;
      char** result;


      rc = sqlite3_get_table(db, zSql, &result, &iRow, &iCol, 0);

      if(rc == SQLITE_OK)
	{
	  curIndex = atoi(result[1]);
	  // curIndex++;
	}
      else
	{
	  cout<<"Can't execute the SQL statement"<<zSql<<endl;

	}

      sprintf(zSql, "select * from msScan, msScanData where startScanNumber=%d "
	      "AND id=scanID", scNum);
      if(!executeSqlStmt(s,zSql))
	cout<<scNum<<" can't be found in the database!"<<endl;

    }
  else
    {
      while(true)
	{
	  curIndex++;
	  if(curIndex > lastIndex)
	    return false;

	  sprintf(zSql, "select * from msScan, msScanData where id=%d "
		  "AND id=scanID",curIndex);
	  if(executeSqlStmt(s,zSql))
	    break;
	}


    }

  return true;

}

bool MSReader::executeSqlStmt(Spectrum& s, char* zSql)
{
  bool isSameLevel=false;
  sqlite3_stmt *pStmt;
  int rc;

  rc = sqlite3_prepare(db, zSql, -1, &pStmt, 0);

  if( rc!=SQLITE_OK ){
    cout<<"can't prepare SQL statement!"<<rc<<endl;
    exit(1);
  }

  rc = sqlite3_step(pStmt);

  char actMethod[1024];
  int charge=-1;
  int msLevel=-1;
  int scanID;
  MSSpectrumType msType;

  if( rc==SQLITE_ROW ){

    //add filter if filter set
    msLevel = sqlite3_column_int(pStmt,4);
    switch(msLevel){
    case 1:
      msType = MS1;
      break;
    case 2:
      msType = MS2;
      break;
    case 3:
      msType = MS3;
      break;
    default:
      break;
    }

    if(find(filter.begin(),filter.end(),msType) != filter.end())
      {
	scanID=sqlite3_column_int(pStmt,0);
	s.setScanID(sqlite3_column_int(pStmt,0));

	s.setScanNumber(sqlite3_column_int(pStmt,2));
	s.setScanNumber(sqlite3_column_int(pStmt,3),true);
	s.setMsLevel(msLevel);
	s.setMZ(sqlite3_column_double(pStmt,5));
	s.setCharge(sqlite3_column_int(pStmt,6));
	readChargeTable(scanID,s);
	s.setRTime((float)sqlite3_column_double(pStmt,9));

	strcpy(actMethod,reinterpret_cast<const char*>(sqlite3_column_text(pStmt,10)));
	if(strcmp(actMethod,"CID") == 0)
	  s.setActivationMethod(CID);
	if(strcmp(actMethod, "ETD") == 0)
	  s.setActivationMethod(ETD);
	if(strcmp(actMethod, "ECD") == 0)
	  s.setActivationMethod(ECD);
	if(strcmp(actMethod, "PQD") == 0)
	  s.setActivationMethod(PQD);
	if(strcmp(actMethod, "HCD") == 0)
	  s.setActivationMethod(HCD);
	if(strcmp(actMethod, "UNKNOWN") == 0)
	  s.setActivationMethod(na);

	int numPeaks = sqlite3_column_int(pStmt,12);

	int numBytes1=sqlite3_column_bytes(pStmt,14);
	unsigned char* comprM = (unsigned char*)sqlite3_column_blob(pStmt,14);
	int numBytes2=sqlite3_column_bytes(pStmt,15);
	unsigned char* comprI = (unsigned char*)sqlite3_column_blob(pStmt,15);


	getUncompressedPeaks(s,numPeaks, numBytes1,comprM, numBytes2,comprI);
	isSameLevel = true;
      }

  }
  rc = sqlite3_finalize(pStmt);
  return isSameLevel;

}

void MSReader::readChargeTable(int scanID, Spectrum& s)
{
  char zSql[8192];
  sprintf(zSql, "select charge, mass from MS2FileScanCharge where scanID=%d", scanID);
  sqlite3_stmt *pStmt;
  int rc;

  rc = sqlite3_prepare(db, zSql, -1, &pStmt, 0);

  if( rc!=SQLITE_OK ){
    cout<<"can't prepare SQL statement!"<<rc<<endl;
    exit(1);
  }

  rc = sqlite3_step(pStmt);
  int charge;
  double MH;
  while(rc == SQLITE_ROW)
    {
      charge = sqlite3_column_int(pStmt,0);
      MH = sqlite3_column_double(pStmt,1);

      s.addZState(charge, MH);
      rc=sqlite3_step(pStmt);
    }
  rc = sqlite3_finalize(pStmt);
}



void MSReader::getUncompressedPeaks(Spectrum& s, int& numPeaks, int& mzLen, unsigned char* comprM, int& intensityLen, unsigned char* comprI)
{
  int i;


  //variables for compressed files
  uLong uncomprLen;
  double *mz;
  float *intensity;

  mz = new double[numPeaks];
  uncomprLen=numPeaks*sizeof(double);
  uncompress((Bytef*)mz, &uncomprLen, comprM, mzLen);

  intensity = new float[numPeaks];
  uncomprLen=numPeaks*sizeof(float);
  uncompress((Bytef*)intensity, &uncomprLen, comprI, intensityLen);

  for(i=0;i<numPeaks;i++){
    s.add(mz[i],intensity[i]);

  }
  delete [] mz;
  delete [] intensity;
}

void MSReader::sql_stmt(const char* stmt)
{

  char *errmsg;
  int   ret;

  ret = sqlite3_exec(db, stmt, 0, 0, &errmsg);

  if (ret != SQLITE_OK)
    {
      printf("Error in statement: %s [%s].\n", stmt, errmsg);
    }
}

#endif

bool MSReader::readFile(const char* c, MSFileFormat f, Spectrum& s, int scNum){

  //Redirect functions to appropriate places, if possible.
  lastFileFormat = f;
  switch(f){
  case ms1:
  case ms2:
  case  zs:
  case uzs:
    return readFile(c,true,s,scNum);
    break;
  case bms1:
  case bms2:
    setCompression(false);
    return readFile(c,false,s,scNum);
    break;
  case cms1:
  case cms2:
    setCompression(true);
    return readFile(c,false,s,scNum);
    break;
  case mzXML:
  case mzData:
    break;
  case raw:
    #ifdef _MSC_VER
    //only read the raw file if the dll was present and loaded.
      if(bRaw) return readRawFile(c,s,scNum);
      else return false;
    #else
      //raw file is not supported
      return false;
    #endif
    break;
  case sqlite:
  case psm:
    #ifndef _NOSQLITE
    return readSqlite(c,s,scNum);
    #else
    //sqlite support disabled
    return false;
    #endif
    break;
  case dunno:
  default:
    return false;
    break;
  }

	//if we got here, it's because we're reading mzXML format

	ramp_fileoffset_t indexOffset;
	ScanHeaderStruct scanHeader;
	RAMPREAL *pPeaks;
	int i,j;

	if(c!=NULL) {
		//open the file if new file was requested
		if(rampFileOpen) {
			rampCloseFile(rampFileIn);
			rampFileOpen=false;
			free(pScanIndex);
		};
		rampFileIn = rampOpenFile(c);
		if (rampFileIn == NULL) {
      cerr << "ERROR: Failure reading input file " << c << endl;
      return false;
		}
		rampFileOpen=true;

		//read the index
		indexOffset = getIndexOffset(rampFileIn);
		pScanIndex = readIndex(rampFileIn,indexOffset,&rampLastScan);
		rampIndex=0;

	} else {
		//if no new file requested, check to see if one is open already
		if (rampFileIn == NULL) return false;
	}


	//clear any spectrum data
	s.clear();

	MSSpectrumType mslevel;

	//read scan header
	if(scNum!=0) {
    rampIndex=scNum;
    /* Henry Lam fixed ramp to take scan numbers as indexes,
       see ramp.cpp commments marked HENRY
		rampIndex=0;
		for(i=1;i<rampLastScan;i++){
			readHeader(rampFileIn, pScanIndex[i], &scanHeader);
			if(scanHeader.acquisitionNum==scNum) {
				rampIndex=i;
				break;
			};
		};
		if(rampIndex==0) return false;
    */
    // Make sure Henry's code works as expected.  Consider making this debug only.
    readHeader(rampFileIn, pScanIndex[rampIndex], &scanHeader);
    if (scanHeader.acquisitionNum != scNum && scanHeader.acquisitionNum != -1)
    {
      cerr << "ERROR: Failure reading scan, index corrupted.  Line endings may have changed during transfer." << flush;
      exit(1);
    }
		/*
		switch(filter){
		case MS1:
			if(scanHeader.msLevel!=1)	return false;
			break;
		case MS2:
			if(scanHeader.msLevel!=2)	return false;
			break;
		case MS3:
			if(scanHeader.msLevel!=3)	return false;
			break;
		default:
			//no filter
			break;
		};
		*/
		switch(scanHeader.msLevel){
		case 1:
		  mslevel = MS1;
		  break;
		case 2:
		  mslevel = MS2;
		  break;
		case 3:
		  mslevel = MS3;
		  break;
		default:
		  break;
		}

		if(find(filter.begin(), filter.end(), mslevel) != filter.end())
		{
		  s.setMsLevel(scanHeader.msLevel);
		  s.setScanNumber(scanHeader.acquisitionNum);
		  s.setScanNumber(scanHeader.acquisitionNum,true);
		  s.setRTime((float)scanHeader.retentionTime);
		  if(strlen(scanHeader.activationMethod)>1){
		    if(strcmp(scanHeader.activationMethod,"CID")==0) s.setActivationMethod(CID);
          else if(strcmp(scanHeader.activationMethod,"ECD")==0) s.setActivationMethod(ECD);
          else if(strcmp(scanHeader.activationMethod,"ETD")==0) s.setActivationMethod(ETD);
          else if(strcmp(scanHeader.activationMethod,"PQD")==0) s.setActivationMethod(PQD);
          else if(strcmp(scanHeader.activationMethod,"HCD")==0) s.setActivationMethod(HCD);
		    else s.setActivationMethod(na);
		  };
		  if(scanHeader.msLevel>1) s.setMZ(scanHeader.precursorMZ);
		  if(scanHeader.precursorCharge>0) s.addZState(scanHeader.precursorCharge,scanHeader.precursorMZ*scanHeader.precursorCharge-(scanHeader.precursorCharge-1)*1.00727649);
		  pPeaks = readPeaks(rampFileIn, pScanIndex[rampIndex]);
		  j=0;
		  for(i=0;i<scanHeader.peaksCount;i++){
		  	s.add((double)pPeaks[j],(float)pPeaks[j+1]);
			  j+=2;
		  }
		}
		else
		  return false;

  } else /* if scnum == 0 */ {

		//read next index
	  while(true){
	    rampIndex++;

	    //reached end of file
	    if(rampIndex>rampLastScan) {
				//rampCloseFile(rampFileIn);
				//rampFileIn = NULL;
				return false;
			};
			readHeader(rampFileIn, pScanIndex[rampIndex], &scanHeader);
			/*
			switch(filter){
			case MS1:
				if(scanHeader.msLevel!=1)	continue;
				break;
			case MS2:
				if(scanHeader.msLevel!=2)	continue;
				break;
			case MS3:
				if(scanHeader.msLevel!=3)	continue;
				break;
			default:
				//no filter
				break;
			};
			*/
			switch(scanHeader.msLevel){
			case 1:
			  mslevel = MS1;
			  break;
			case 2:
			  mslevel = MS2;
			  break;
			case 3:
			  mslevel = MS3;
			  break;
			default:
			  break;
			}
			if(find(filter.begin(), filter.end(), mslevel) != filter.end())
			  break;
			//if we got here, we passed the filter.

		};
		s.setMsLevel(scanHeader.msLevel);
		s.setScanNumber(scanHeader.acquisitionNum);
		s.setScanNumber(scanHeader.acquisitionNum,true);
		s.setRTime((float)scanHeader.retentionTime);
		if(strlen(scanHeader.activationMethod)>1){
		  if(strcmp(scanHeader.activationMethod,"CID")==0) s.setActivationMethod(CID);
        else if(strcmp(scanHeader.activationMethod,"ECD")==0) s.setActivationMethod(ECD);
        else if(strcmp(scanHeader.activationMethod,"ETD")==0) s.setActivationMethod(ETD);
        else if(strcmp(scanHeader.activationMethod,"PQD")==0) s.setActivationMethod(PQD);
        else if(strcmp(scanHeader.activationMethod,"HCD")==0) s.setActivationMethod(HCD);
		  else s.setActivationMethod(na);
		};
		if(scanHeader.msLevel>1) s.setMZ(scanHeader.precursorMZ);
		if(scanHeader.precursorCharge>0) s.addZState(scanHeader.precursorCharge,scanHeader.precursorMZ*scanHeader.precursorCharge-(scanHeader.precursorCharge-1)*1.00727649);
		pPeaks = readPeaks(rampFileIn, pScanIndex[rampIndex]);
		j=0;
		for(i=0;i<scanHeader.peaksCount;i++){
			s.add((double)pPeaks[j],(float)pPeaks[j+1]);
			j+=2;
		};

    }

	free(pPeaks);
	return true;

};

void MSReader::setFilter(MSSpectrumType m){
  filter.clear();
  filter.push_back(m);
}

void MSReader::setFilter(vector<MSSpectrumType>& m){
  for(unsigned int i=0; i<m.size();i++)
    filter.push_back(m.at(i));
}

void MSReader::setCompression(bool b){
	compressMe=b;
}


void MSReader::setAverageRaw(bool b, int width, long cutoff){
  #ifdef _MSC_VER
  rawAvg=b;
  rawAvgWidth=width;
  rawAvgCutoff=cutoff;
  #endif
}

void MSReader::setLabel(bool b){
  #ifdef _MSC_VER
  rawLabel=b;
  #endif
}

void MSReader::setRawFilter(char *c){
  #ifdef _MSC_VER
  strcpy(rawUserFilter,c);
  #endif
}

void MSReader::setHighResMGF(bool b){
  highResMGF=b;
}

void MSReader::writeCompressSpec(FILE* fileOut, Spectrum& s){

	int j;

	//file compression
	int err;
	uLong len;
	unsigned char *comprM, *comprI;
  uLong comprLenM, comprLenI;
	double *pD;
	float *pF;
	uLong sizeM;
	uLong sizeI;

	//Build arrays to hold scan prior to compression
	// Ideally, we would just use the scan vectors, but I don't know how yet.
	pD = new double[s.size()];
	pF = new float[s.size()];
	for(j=0;j<s.size();j++){
		pD[j]=s.at(j).mz;
		pF[j]=s.at(j).intensity;
	};

	//compress mz
	len = (uLong)s.size()*sizeof(double);
	sizeM = len;
	comprLenM = compressBound(len);
	comprM = (unsigned char*)calloc((uInt)comprLenM, 1);
	err = compress(comprM, &comprLenM, (const Bytef*)pD, len);

	//compress intensity
	len = (uLong)s.size()*sizeof(float);
	sizeI = len;
	comprLenI = compressBound(len);
	comprI = (unsigned char*)calloc((uInt)comprLenI, 1);
	err = compress(comprI, &comprLenI, (const Bytef*)pF, len);

	j=(int)comprLenM;
	fwrite(&j,4,1,fileOut);
	j=(int)comprLenI;
	fwrite(&j,4,1,fileOut);
	fwrite(comprM,comprLenM,1,fileOut);
	fwrite(comprI,comprLenI,1,fileOut);

	//clean up memory
	free(comprM);
	free(comprI);
	delete [] pD;
	delete [] pF;

};

void MSReader::readCompressSpec(FILE* fileIn, MSScanInfo& ms, Spectrum& s){

	int i;
	Peak_T p;

	//variables for compressed files
	uLong uncomprLen;
	uLong mzLen, intensityLen;
	unsigned char *compr;
	double *mz;
	float *intensity;

	fread(&i,4,1,fileIn);
	mzLen = (uLong)i;
	fread(&i,4,1,fileIn);
	intensityLen = (uLong)i;

	compr = new unsigned char[mzLen];
	mz = new double[ms.numDataPoints];
	uncomprLen=ms.numDataPoints*sizeof(double);
	fread(compr,mzLen,1,fileIn);
	uncompress((Bytef*)mz, &uncomprLen, compr, mzLen);
	delete [] compr;

	compr = new unsigned char[intensityLen];
	intensity = new float[ms.numDataPoints];
	uncomprLen=ms.numDataPoints*sizeof(float);
	fread(compr,intensityLen,1,fileIn);
	uncompress((Bytef*)intensity, &uncomprLen, compr, intensityLen);
	delete [] compr;

	for(i=0;i<ms.numDataPoints;i++){
		p.mz = mz[i];
		p.intensity = intensity[i];
		s.add(p);
	};

	delete [] mz;
	delete [] intensity;

};

void MSReader::writeTextSpec(FILE* fileOut, Spectrum& s) {

	int i,j,k;
	char t[64];

  if(exportMGF){
    //MGF spectrum header is here
    if(highResMGF){
      for(i=0;i<s.sizeZ();i++){
        fprintf(fileOut,"BEGIN IONS\n");
        fprintf(fileOut,"PEPMASS=%.*f\n",6,s.atZ(i).mz);
        fprintf(fileOut,"CHARGE=%d+\n",s.atZ(i).z);
        fprintf(fileOut,"TITLE=%s.%d.%d.%d %d %.4f\n","test",s.getScanNumber(),s.getScanNumber(true),s.atZ(i).z,i,s.getRTime());
        for(j=0;j<s.size();j++){
		      sprintf(t,"%.*f",iIntensityPrecision,s.at(j).intensity);
		      k=strlen(t);
		      if(k>2 && iIntensityPrecision>0){
		        if(t[0]=='0'){
		          fprintf(fileOut,"%.*f 0\n",iMZPrecision,s.at(j).mz);
			      } else if(t[k-1]=='0'){
			        fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision-1,s.at(j).intensity);
       			} else {
			        fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
			      }
		      } else {
			      fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
		      }
	      }
        fprintf(fileOut,"END IONS\n");
      }

    } else {
      fprintf(fileOut,"BEGIN IONS\n");
      fprintf(fileOut,"PEPMASS=%.*f\n",6,s.getMZ());
      if(s.sizeZ()==1){
        if(s.atZ(0).z==1) fprintf(fileOut,"CHARGE=1+\n");
      }
      fprintf(fileOut,"TITLE=%s.%d.%d.%d %d %.4f\n","test",s.getScanNumber(),s.getScanNumber(true),s.atZ(0).z,0,s.getRTime());
      for(j=0;j<s.size();j++){
		    sprintf(t,"%.*f",iIntensityPrecision,s.at(j).intensity);
		    k=strlen(t);
		    if(k>2 && iIntensityPrecision>0){
		      if(t[0]=='0'){
		        fprintf(fileOut,"%.*f 0\n",iMZPrecision,s.at(j).mz);
			    } else if(t[k-1]=='0'){
			      fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision-1,s.at(j).intensity);
      		} else {
			      fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
			    }
		    } else {
			    fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
		    }
	    }
      fprintf(fileOut,"END IONS\n");
    }
    return;
  }

  //Only use this code if not writing MGF file
	for(j=0;j<s.size();j++){
		sprintf(t,"%.*f",iIntensityPrecision,s.at(j).intensity);
		k=strlen(t);
		if(k>2 && iIntensityPrecision>0){
			if(t[0]=='0'){
				fprintf(fileOut,"%.*f 0\n",iMZPrecision,s.at(j).mz);
			} else if(t[k-1]=='0'){
				fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision-1,s.at(j).intensity);
			} else {
				fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
			}
		} else {
			fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
		}
	}

}

void MSReader::writeBinarySpec(FILE* fileOut, Spectrum& s) {
	int j;

	for(j=0;j<s.size();j++){
		fwrite(&s.at(j).mz,8,1,fileOut);
		fwrite(&s.at(j).intensity,4,1,fileOut);
	};

};

void MSReader::writeSpecHeader(FILE* fileOut, bool text, Spectrum& s) {

	//MSScanInfo ms;
  double d;
  float f;
  int i;
	MSSpectrumType mft;
	int j;

	//output scan info
	if(text){

    //MSx spectrum header is here
		mft=s.getFileType();
   	if(mft==MS2 || mft==MS3 || mft==SRM){
    	fprintf(fileOut,"S\t%d\t%d\t%.*f\n",s.getScanNumber(),s.getScanNumber(true),4,s.getMZ());
		} else {
	  	fprintf(fileOut,"S\t%d\t%d\n",s.getScanNumber(),s.getScanNumber(true));
		}
  	if(s.getRTime()>0) fprintf(fileOut,"I\tRTime\t%.*f\n",4,s.getRTime());
    if(s.getBPI()>0) fprintf(fileOut,"I\tBPI\t%.*f\n",2,s.getBPI());
    if(s.getBPM()>0) fprintf(fileOut,"I\tBPM\t%.*f\n",4,s.getBPM());
    if(s.getConversionA()!=0) fprintf(fileOut,"I\tConvA\t%.*f\n",4,s.getConversionA());
    if(s.getConversionB()!=0) fprintf(fileOut,"I\tConvB\t%.*f\n",4,s.getConversionB());
    if(s.getTIC()>0) fprintf(fileOut,"I\tTIC\t%.*f\n",2,s.getTIC());
    if(s.getIonInjectionTime()>0) fprintf(fileOut,"I\tIIT\t%.*f\n",4,s.getIonInjectionTime());
    for(j=0;j<s.sizeEZ();j++){
      fprintf(fileOut,"I\tEZ\t%d\t%.*f\t%.*f\t%.*f\n",s.atEZ(j).z,4,s.atEZ(j).mh,4,s.atEZ(j).pRTime,1,s.atEZ(j).pArea);
  	}
	  for(j=0;j<s.sizeZ();j++){
		 	fprintf(fileOut,"Z\t%d\t%.*f\n",s.atZ(j).z,4,s.atZ(j).mz);
		}

	} else {
    i=s.getScanNumber();
    fwrite(&i,4,1,fileOut);

    i=s.getScanNumber(true);
    fwrite(&i,4,1,fileOut);

    d=s.getMZ();
    fwrite(&d,8,1,fileOut);

    f=s.getRTime();
    fwrite(&f,4,1,fileOut);

    f=s.getBPI();
    fwrite(&f,4,1,fileOut);

    d=s.getBPM();
    fwrite(&d,8,1,fileOut);

    d=s.getConversionA();
    fwrite(&d,8,1,fileOut);

    d=s.getConversionB();
    fwrite(&d,8,1,fileOut);

    d=s.getTIC();
    fwrite(&d,8,1,fileOut);

    f=s.getIonInjectionTime();
    fwrite(&f,4,1,fileOut);

    i=s.sizeZ();
    fwrite(&i,4,1,fileOut);

    i=s.sizeEZ();
    fwrite(&i,4,1,fileOut);

    i=s.size();
    fwrite(&i,4,1,fileOut);
    /*
		ms.scanNumber[0]=ms.scanNumber[1]=s.getScanNumber();
		ms.rTime=s.getRTime();
		ms.numDataPoints=s.size();
		ms.numZStates=s.sizeZ();
		fwrite(&ms,sizeof(MSScanInfo),1,fileOut);
    */
    for(j=0;j<s.sizeZ();j++){
			fwrite(&s.atZ(j).z,4,1,fileOut);
			fwrite(&s.atZ(j).mz,8,1,fileOut);
		}

    for(j=0;j<s.sizeEZ();j++){
			fwrite(&s.atEZ(j).z,4,1,fileOut);
			fwrite(&s.atEZ(j).mh,8,1,fileOut);
      fwrite(&s.atEZ(j).pRTime,4,1,fileOut);
      fwrite(&s.atEZ(j).pArea,4,1,fileOut);
		}
	}

}

void MSReader::readSpecHeader(FILE *fileIn, MSScanInfo &ms){

  fread(&ms.scanNumber[0],4,1,fileIn);
  if(feof(fileIn)) return;
  fread(&ms.scanNumber[1],4,1,fileIn);
  fread(&ms.mz,8,1,fileIn);
  fread(&ms.rTime,4,1,fileIn);

  if(iVersion>=2){
    fread(&ms.BPI,4,1,fileIn);
    fread(&ms.BPM,8,1,fileIn);
    fread(&ms.convA,8,1,fileIn);
    fread(&ms.convB,8,1,fileIn);
    fread(&ms.TIC,8,1,fileIn);
    fread(&ms.IIT,4,1,fileIn);
  }

  fread(&ms.numZStates,4,1,fileIn);

  if(iVersion>=3) fread(&ms.numEZStates,4,1,fileIn);
  else ms.numEZStates=0;

  fread(&ms.numDataPoints,4,1,fileIn);

}

MSFileFormat MSReader::checkFileFormat(const char *fn){

  int i;

  //check extension first - we must trust MS1 & MS2 & ZS & UZS
  i=strlen(fn);
  if(strcmp(fn+(i-4),".ms1")==0 || strcmp(fn+(i-4),".MS1")==0 ) return ms1;
  if(strcmp(fn+(i-4),".ms2")==0 || strcmp(fn+(i-4),".MS2")==0 ) return ms2;
  if(strcmp(fn+(i-3),".zs")==0 || strcmp(fn+(i-3),".ZS")==0 ) return zs;
  if(strcmp(fn+(i-4),".uzs")==0 || strcmp(fn+(i-4),".UZS")==0 ) return uzs;
  if(strcmp(fn+(i-6),".msmat")==0 || strcmp(fn+(i-6),".MSMAT")==0 ) return msmat_ff;

  //RAW could mean anything, Thermo should know better
  if(strcmp(fn+(i-4),".raw")==0 || strcmp(fn+(i-4),".RAW")==0 ) {
    #ifdef _MSC_VER
    if(bRaw) {
      return raw;
    } else {
      cout << "Thermo RAW files not suppored on this computer. Please check for Xcalibur installation." << endl;
      return dunno;
    }
    #else
      cout << "Thermo RAW files are only supported on Windows operating system with Xcalibur installed." << endl;
      return dunno;
    #endif
  }

  //For now, trust mzXML & mzData also
  if(strcmp(fn+(i-6),".mzXML")==0 || strcmp(fn+(i-6),".mzxml")==0 || strcmp(fn+(i-6),".MZXML")==0 ) return mzXML;
  if(strcmp(fn+(i-7),".mzData")==0 || strcmp(fn+(i-7),".mzdata")==0 || strcmp(fn+(i-7),".MZDATA")==0 ) return mzData;

  //MGF format
  if(strcmp(fn+(i-4),".mgf")==0 || strcmp(fn+(i-4),".MGF")==0 ) return mgf;

  //add the sqlite3 format
  if(strcmp(fn+(i-8),".sqlite3")==0 || strcmp(fn+(i-8),".SQlite3")==0 || strcmp(fn+(i-8),".SQLite3")==0 ) return sqlite;

  if(strcmp(fn+(i-4), ".psm") == 0 || strcmp(fn+(i-4), ".PSM") == 0) return psm;

  //We can check headers for other formats
  FILE *f;
  f=fopen(fn,"rb");
  if(f==NULL) {
    cout << "Error! Cannot open file to check file format." << endl;
    return dunno;
  }
  fread(&i,4,1,f);
  fread(&iVersion,4,1,f);
  fclose(f);

  //cout << fn << "  FType = " << i << endl;
  switch(i){
    case 1: return bms1;
    case 2: return cms1;
    case 3: return bms2;
    case 4: return cms2;
    default: return dunno;
  }

  return dunno;

}

#ifdef _MSC_VER
bool MSReader::readRawFile(const char *c, Spectrum &s, int scNum){
  HRESULT lRet;
  TCHAR pth[MAX_PATH];
  VARIANT varMassList;
	VARIANT varPeakFlags;
  SAFEARRAY FAR* psa;
  DataPeak* pDataPeaks = NULL;
	long lArraySize=0;
	double dRTime;
	double precursorMZ,highmass=0.0;
	MSSpectrumType MSn;
	int Charge;
  int rawCharge;
  double rawMZ;
  long j;
  long i;
  int k;
  double pw;
  double pm1;
  double d;
  char chFilter[256];
  char curFilter[256];
  strcpy(chFilter,"");
  strcpy(curFilter,"");

  //For gathering averaged scans
  long FirstBkg1=0;
  long LastBkg1=0;
  long FirstBkg2=0;
  long LastBkg2=0;
  int widthCount;
  long lowerBound;
  long upperBound;
  BSTR rawFilter=NULL;

  //Additional Scan Information
  VARIANT ConversionA;
  VARIANT ConversionB;
  double TIC;
  VARIANT IIT;  //ion injection time
  double BPM;   //Base peak mass
  double BPI;   //Base peak intensity
  long tl;      //temp long value
  double td;    //temp double value

  bool bCheckNext;

  double* pdval;

  s.clear();

  if(c==NULL){
    if(scNum>0) rawCurSpec=scNum;
    else rawCurSpec++;
    if(rawCurSpec>rawTotSpec) {
      if(rawFileOpen) lRet = m_Raw->Close();
      rawFileOpen=false;
      return false;
    }
  } else {
    if(rawFileOpen) {
      lRet = m_Raw->Close();
      rawFileOpen=false;
    }
    MultiByteToWideChar(CP_ACP,0,c,-1,(LPWSTR)pth,MAX_PATH);
    lRet = m_Raw->Open((LPWSTR)pth);
	  if(lRet != ERROR_SUCCESS) return false;
	  else lRet = m_Raw->SetCurrentController(0,1);
    rawFileOpen=true;
    m_Raw->GetNumSpectra(&rawTotSpec);
    if(scNum>0) rawCurSpec=scNum;
    else rawCurSpec=1;
    if(rawCurSpec>rawTotSpec) {
      if(rawFileOpen) lRet = m_Raw->Close();
      rawFileOpen=false;
      return false;
    }
  }

	//Initialize structures for excalibur
	VariantInit(&varMassList);
	VariantInit(&varPeakFlags);

  //if the filter was set, make sure we pass the filter
  while(true){
    //cout << rawCurSpec << endl;
	  MSn = EvaluateFilter(rawCurSpec, &precursorMZ, curFilter, rawCharge, rawMZ);

    //check for spectrum filter
    if(strlen(rawUserFilter)>0){
      bCheckNext=false;
      if(rawUserFilterExact) {
        if(strcmp(curFilter,rawUserFilter)!=0)bCheckNext=true;
      } else {
        if(strstr(curFilter,rawUserFilter)==NULL) bCheckNext=true;
      }

      //if string doesn't match, get next scan until it does match or EOF
      if(bCheckNext){
        if(scNum>0) return false;
        rawCurSpec++;
        if(rawCurSpec>rawTotSpec) {
          if(rawFileOpen) lRet = m_Raw->Close();
          rawFileOpen=false;
          return false;
        }
        continue;
      }
    }

    //check for charge state filter
    if(filter.size()>0 && find(filter.begin(), filter.end(), MSn) == filter.end()) {
      if(scNum>0) return false;
      rawCurSpec++;
      if(rawCurSpec>rawTotSpec) {
        if(rawFileOpen) lRet = m_Raw->Close();
        rawFileOpen=false;
        return false;
      }
    } else {
      break;
    }
  }

  //Get any additional information
  VariantInit(&IIT);
  VariantInit(&ConversionA);
  VariantInit(&ConversionB);
  m_Raw->GetTrailerExtraValueForScanNum(rawCurSpec, _T("Ion Injection Time (ms):") , &IIT);
  m_Raw->GetTrailerExtraValueForScanNum(rawCurSpec, _T("Conversion Parameter A:") , &ConversionA);
  m_Raw->GetTrailerExtraValueForScanNum(rawCurSpec, _T("Conversion Parameter B:") , &ConversionB);
  m_Raw->GetScanHeaderInfoForScanNum(rawCurSpec, &tl, &td, &td, &td, &TIC, &BPM, &BPI, &tl, &tl, &td);
  m_Raw->RTFromScanNum(rawCurSpec,&dRTime);

  //Get the peaks
  if(rawAvg){
    widthCount=0;
    lowerBound=0;
    upperBound=0;
    for(i=rawCurSpec-1;i>0;i--){
      EvaluateFilter(i, &d, chFilter, k, d);
      if(strcmp(curFilter,chFilter)==0){
        widthCount++;
        if(widthCount==rawAvgWidth) {
          lowerBound=i;
          break;
        }
      }
    }
    if(lowerBound==0) lowerBound=rawCurSpec; //this will have "edge" effects

    widthCount=0;
    for(i=rawCurSpec+1;i<rawTotSpec;i++){
      EvaluateFilter(i, &d, chFilter, k, d);
      if(strcmp(curFilter,chFilter)==0){
        widthCount++;
        if(widthCount==rawAvgWidth) {
          upperBound=i;
          break;
        }
      }
    }
    if(upperBound==0) upperBound=rawCurSpec; //this will have "edge" effects

    m_Raw->GetFilterForScanNum(i, &rawFilter);
    //cout << lowerBound << " xx " << upperBound << endl;
    j=m_Raw->GetAverageMassList(&lowerBound, &upperBound, &FirstBkg1, &LastBkg1, &FirstBkg2, &LastBkg2,
      rawFilter, 1, rawAvgCutoff, 0, FALSE, &pw, &varMassList, &varPeakFlags, &lArraySize );
    SysFreeString(rawFilter);
    rawFilter=NULL;

  } else {
    if(rawLabel) {
      j=m_Raw->GetLabelData(&varMassList, &varPeakFlags, &rawCurSpec);
    } else {
      j=m_Raw->GetMassListFromScanNum(&rawCurSpec,_T(""),0,0,0,FALSE,&pw,&varMassList,&varPeakFlags,&lArraySize);
    }
  }

	//Handle MS2 and MS3 files differently to create Z-lines
	if(MSn==MS2 || MSn==MS3){

    if(rawCharge>0){
      if(rawMZ>0.0) pm1 = CalcPepMass(rawCharge, rawMZ);
      else pm1 = CalcPepMass(rawCharge, precursorMZ);
      s.addZState(rawCharge,pm1);
    } else {
      Charge = CalcChargeState(precursorMZ, highmass, &varMassList, lArraySize);

      //Charge greater than 0 means the charge state is known
      if(Charge>0){
        pm1 = CalcPepMass(Charge, precursorMZ);
  	    s.addZState(Charge,pm1);

      //Charge of 0 means unknown charge state, therefore, compute +2 and +3 states.
      } else {
        pm1 = CalcPepMass(2, precursorMZ);
        s.addZState(2,pm1);
        pm1 = CalcPepMass(3, precursorMZ);
        s.addZState(3,pm1);
      }

    }

  } //endif MS2 and MS3

	//Set basic scan info
  s.setRawFilter(curFilter);
	s.setScanNumber((int)rawCurSpec);
  s.setScanNumber((int)rawCurSpec,true);
	s.setRTime((float)dRTime);
	s.setFileType(MSn);
  s.setBPI((float)BPI);
  s.setBPM(BPM);
  s.setConversionA(ConversionA.dblVal);
  s.setConversionB(ConversionB.dblVal);
  s.setTIC(TIC);
  s.setIonInjectionTime(IIT.fltVal);
	if(MSn==MS2 || MSn==MS3 || MSn==SRM) s.setMZ(precursorMZ);
  switch(MSn){
    case MS1: s.setMsLevel(1); break;
    case MS2: s.setMsLevel(2); break;
    case MS3: s.setMsLevel(3); break;
    default: s.setMsLevel(0); break;
  }

  //Clean up memory
  VariantClear(&IIT);
  VariantClear(&ConversionA);
  VariantClear(&ConversionB);

  if(rawLabel){
    psa = varMassList.parray;
    lArraySize = psa->rgsabound[0].cElements;
    pdval = (double *) psa->pvData;
    for(j=0;j<lArraySize;j++) s.add((double)pdval[j*6],(float)pdval[j*6+1]);
  } else {
	  psa = varMassList.parray;
	  SafeArrayAccessData( psa, (void**)(&pDataPeaks) );
	  for(j=0;j<lArraySize;j++)	s.add(pDataPeaks[j].dMass,(float)pDataPeaks[j].dIntensity);
	  SafeArrayUnaccessData( psa );
  }

	//Clear Xcalibur structures
	VariantClear(&varMassList);
	VariantClear(&varPeakFlags);

  return true;

}

MSSpectrumType MSReader::EvaluateFilter(long scan, double *precursormz, char* chFilter, int &thermoCharge, double &thermoMZ) {

  USES_CONVERSION;

  BSTR Filter = NULL;
	string cFilter;
	string tStr;
	int start,stop;
  bool bSRM;
  int i;
  VARIANT RAWCharge;
  VARIANT RAWmz;

  //Initialize raw values to default
  thermoCharge=0;
  thermoMZ=0.0;

	m_Raw->GetFilterForScanNum(scan, &Filter);
  strcpy(chFilter,W2A(Filter));
  cFilter=chFilter;

  //Check for SRM tag - this isn't used yet, but keeping it for future needs
  bSRM=false;
  stop=cFilter.find("SRM");
  if(stop>-1)bSRM=true;

	//search for ms2/ms3/msn
	stop=cFilter.find("@");
  if(stop > -1){
    //GF:moved SRM ms2 to the else block, as it does not have a '@' symbol (with data generated 05/18/08 on quantum)
    //MH:Data generated 09/06/07 on the quantum has both SRM and '@' symbol. Will check for SRM in both blocks

		start=cFilter.find("Full ms2",0);
		if(start>-1){
			start+=9;
			tStr=cFilter.substr(start,stop-start);
			*precursormz=atof(&tStr[0]);
      VariantInit(&RAWCharge);
      VariantInit(&RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Monoisotopic M/Z:") , &RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Charge State:") , &RAWCharge);
      thermoCharge=RAWCharge.iVal;
      thermoMZ=RAWmz.dblVal;
      VariantClear(&RAWmz);
      VariantClear(&RAWCharge);
			return MS2;
    }

    start=cFilter.find("SRM ms2",0);
		if(start>-1){
			start+=9;
			tStr=cFilter.substr(start,stop-start);
			*precursormz=atof(&tStr[0]);
			return SRM;
		} else {
			start=cFilter.find(" ",stop+1);
			stop=cFilter.find("@",stop+1);
			tStr=cFilter.substr(start,stop-start);
			*precursormz=atof(&tStr[0]);
			return MS3;
	  }

	//search for Ultra/Zoom scans,
  //and SRM scans (GF 05/18/08)
  } else {
		if(cFilter.find("u Z ms")>-1) {
      VariantInit(&RAWCharge);
      VariantInit(&RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Monoisotopic M/Z:") , &RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Charge State:") , &RAWCharge);
      thermoCharge=RAWCharge.iVal;
      thermoMZ=RAWmz.dblVal;
      VariantClear(&RAWmz);
      VariantClear(&RAWCharge);
			return UZS;
    }
		if(cFilter.find("Z ms")>-1){
      VariantInit(&RAWCharge);
      VariantInit(&RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Monoisotopic M/Z:") , &RAWmz);
      m_Raw->GetTrailerExtraValueForScanNum(scan, _T("Charge State:") , &RAWCharge);
      thermoCharge=RAWCharge.iVal;
      thermoMZ=RAWmz.dblVal;
      VariantClear(&RAWmz);
      VariantClear(&RAWCharge);
			return ZS;
    }

    i=cFilter.find("Full ms");
    if( i > -1)	return MS1;
   
    if(cFilter.find("SRM ms2")>-1) {
			start = cFilter.find("SRM ms2",0);
      stop = cFilter.find("[",start);
		  start+=8;
			tStr=cFilter.substr(start,stop-start);
		  *precursormz=atof(&tStr[0]);
			return SRM;
    }
	}

	//When all else fails, just set Unknown
	return Unspecified;
}

int MSReader::CalcChargeState(double precursormz, double highmass, VARIANT* varMassList, long nArraySize) {
// Assumes spectrum is +1 or +2.  Figures out charge by
// seeing if signal is present above the parent mass
// indicating +2 (by taking ratio above/below precursor)

	bool bFound;
	long i, iStart;
	double dLeftSum,dRightSum;
	double FractionWindow;
	double CorrectionFactor;

	dLeftSum = 0.00001;
	dRightSum = 0.00001;

	DataPeak* pDataPeaks = NULL;
	SAFEARRAY FAR* psa = varMassList->parray;
	SafeArrayAccessData( psa, (void**)(&pDataPeaks) );

//-------------
// calc charge
//-------------
	bFound=false;
	i=0;
	while(i<nArraySize && !bFound){
    if(pDataPeaks[i].dMass < precursormz - 20){
			//do nothing
		} else {
			bFound = true;
      iStart = i;
    }
    i++;
	}
	if(!bFound) iStart = nArraySize;

	for(i=0;i<iStart;i++)	dLeftSum = dLeftSum + pDataPeaks[i].dIntensity;

	bFound=false;
	i=0;
	while(i<nArraySize && !bFound){
    if(pDataPeaks[i].dMass < precursormz + 20){
			//do nothing
		} else {
      bFound = true;
      iStart = i;
    }
    i++;
	}

	if(!bFound) {
		SafeArrayUnaccessData( psa );
		psa = NULL;
		pDataPeaks = NULL;
		return 1;
	}
	if(iStart = 0) iStart++;

	for(i=iStart;i<nArraySize;i++) dRightSum = dRightSum + pDataPeaks[i].dIntensity;

	if(precursormz * 2 < highmass){
    CorrectionFactor = 1;
	} else {
    FractionWindow = (precursormz * 2) - highmass;
    CorrectionFactor = (precursormz - FractionWindow) / precursormz;
	}

	if(dLeftSum > 0 && (dRightSum / dLeftSum) < (0.2 * CorrectionFactor)){
		SafeArrayUnaccessData( psa );
		psa=NULL;
		pDataPeaks=NULL;
		return 1;
	} else {
		SafeArrayUnaccessData( psa );
		psa=NULL;
		pDataPeaks=NULL;
    return 0;  //Set charge to 0 to indicate that both +2 and +3 spectra should be created
	}

  //When all else fails, return 0
  return 0;
}

double MSReader::CalcPepMass(int chargestate, double precursormz){
  double MplusH=0.0;
  MplusH = precursormz * chargestate - ((chargestate-1)*1.00727649);
	return MplusH;
}

void MSReader::setRawFilterExact(bool b){
  rawUserFilterExact=b;
}

bool MSReader::lookupRT(char* c, int scanNum, float& rt){
  HRESULT lRet;
  TCHAR pth[MAX_PATH];
  double drt;

  rt=0.0f;

  if(c!=NULL){
    if(rawFileOpen) {
      lRet = m_Raw->Close();
      rawFileOpen=false;
    }

    MultiByteToWideChar(CP_ACP,0,c,-1,(LPWSTR)pth,MAX_PATH);
    lRet = m_Raw->Open((LPWSTR)pth);
    if(lRet != ERROR_SUCCESS) return false;
    else lRet = m_Raw->SetCurrentController(0,1);
    rawFileOpen=true;
    m_Raw->GetNumSpectra(&rawTotSpec);
    if(scanNum>rawTotSpec) {
      if(rawFileOpen) lRet = m_Raw->Close();
      rawFileOpen=false;
      return false;
    }
  }

  if(m_Raw->RTFromScanNum(scanNum,&drt)==0) {
    rt=(float)drt;
    return true;
  }

  return false;

}

#endif