#include "MSReader.h"
#include <iostream>
using namespace std;

MSReader::MSReader(){
	fileIn=NULL;
	iIntensityPrecision=1;
	iMZPrecision=4;
	filter=Unspecified;
	rampFileOpen=false;
	compressMe=false;
  iFType=0;
  iVersion=0;
};

MSReader::~MSReader(){
	closeFile();
	if(rampFileOpen) {
		rampCloseFile(rampFileIn);
		free(pScanIndex);
	};

};

void MSReader::closeFile(){
  if(fileIn!=NULL) fclose(fileIn);
};

MSHeader& MSReader::getHeader(){
  return header;
};

/* 0 = File opened correctly
   1 = Could not open file
*/
int MSReader::openFile(char *c,bool text){
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
};

/*
Spectrum MSReader::readBinaryFile(char *c, Spectrum& s, int scNum){
	MS1ScanInfo ms;
	Peak_T p;
	int i;

	s.clear();

	if(c!=NULL){
		closeFile();
		if(openFile(c,false)==1) return s;
	} else if(fileIn==NULL) {
		return s;
	};


	fread(&ms,sizeof(MS1ScanInfo),1,fileIn);
	if(scNum!=0){
		while(ms.scanNumber[0]!=scNum){
			fseek(fileIn,ms.numDataPoints*12,1);
			fread(&ms,sizeof(MS1ScanInfo),1,fileIn);
			if(feof(fileIn)) return s;
		};
	};
	if(feof(fileIn)) return s;

	s.setScanNumber(ms.scanNumber[0]);
	s.setRTime(ms.rTime);
	for(i=0;i<ms.numDataPoints;i++){
		fread(&p.mz,8,1,fileIn);
		fread(&p.intensity,4,1,fileIn);
		s.add(p);
	};

	return s;

};
*/

bool MSReader::readFile(char *c, bool text, Spectrum& s, int scNum){
	MSScanInfo ms;
	Peak_T p;
	ZState z;
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
	} else if(fileIn==NULL) {
		return false;
	};

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

				if(compressMe){
					fread(&i,4,1,fileIn);
					mzLen = (uLong)i;
					fread(&i,4,1,fileIn);
					intensityLen = (uLong)i;
					fseek(fileIn,mzLen+intensityLen,1);
				} else {
					fseek(fileIn,ms.numDataPoints*12,1);
				};

				//fread(&ms,sizeof(MSScanInfo),1,fileIn);
        readSpecHeader(fileIn,ms);
				if(feof(fileIn)) return false;
			};
		};
		if(feof(fileIn)) return false;

		//read any charge states (for MS2 files)
		for(i=0;i<ms.numZStates;i++){
			fread(&z.z,4,1,fileIn);
			fread(&z.mz,8,1,fileIn);
			s.addZState(z);
		};

		s.setScanNumber(ms.scanNumber[0]);
    s.setScanNumber(ms.scanNumber[1],true);
		s.setRTime(ms.rTime);
    s.setMZ(ms.mz);

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
			};
		};

		//return success
		return true;

	} else {

		//if reading text files, some parsing is required.
		while(true){
			if(feof(fileIn)) break;
    
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
				};
				break;

			case 'I':
				//I lines are recorded only if they contain retention times
				fgets(tstr,256,fileIn);
				tok=strtok(tstr," \t\n\r");
				tok=strtok(NULL," \t\n\r");
				if(strcmp(tok,"RTime")==0) {
					tok=strtok(NULL," \t\n\r,");
					s.setRTime(atof(tok));
				};
				break;

			case 'S':
				//Scan numbers are recorded and mark all following data is spectrum data
				//until the next tag

				//Reaching an S tag also indicates there are no more header lines
				bDoneHeader=true;

				if(firstScan) {
					//if we are here, a scan was read and we just reached the next scan tag
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
							s.setScanNumber(0);
							if(bScan==false) return false;
							break;
						};
					};
					firstScan=true;
				};
				break;

			case 'Z':
				//Z lines are recorded for MS2 files
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
					};
				};
				//otherwise, read in the line
				fscanf(fileIn,"%lf %f\n",&p.mz,&p.intensity);
				s.add(p);
				break;
			
			default:
				//if the character is not recognized, ignore the entire line.
				fscanf(fileIn,"%s\n",tstr);
				break;
			};
		};

	};
  
  return true;

};


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
  };
  
	fseek(fileIn,lPivot,0);

  return (lFWidth>0 && lPivot>0 && lPivot<lEnd);

};


int MSReader::getPercent(){
  if(fileIn!=NULL){
		return (int)((double)ftell(fileIn)/lEnd*100);
  };
	if(rampFileIn!=NULL){
		return (int)((double)rampIndex/rampLastScan*100);
	};
  return -1;
};

void MSReader::writeFile(char* c, bool text, MSObject& m){

  FILE* fileOut;
  int i;

  //if a filename isn't specified, check to see if the
  //MSObject has a filename.
  if(c == NULL) {
		return;
  } else {
		if(text) fileOut=fopen(c,"wt");
		else fileOut=fopen(c,"wb");
	};

  //output file header lines;
	if(text){
		for(i=0;i<16;i++){
			if(m.getHeader().header[i][0]!='\0') {
				fputs("H\t",fileOut);
				fputs(m.getHeader().header[i],fileOut);
			};
		};
	} else {
    fwrite(&iFType,4,1,fileOut); //file type
    i=1;
    fwrite(&i,4,1,fileOut); //version number - in case we change formats
		fwrite(&m.getHeader(),sizeof(MSHeader),1,fileOut);
	};
  
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
		};

  };
    
	fclose(fileOut);
};

void MSReader::writeFile(char* c, MSFileFormat ff, MSObject& m){

  int i;

  switch(ff){
    case ms1:
    case bms1:
    case cms1:
      for(i=0;i<m.size();i++) m.at(i).setFileType(MS1);
      break;
    case ms2:
    case bms2:
    case cms2:
      for(i=0;i<m.size();i++) m.at(i).setFileType(MS2);
      break;
    case zs:
      for(i=0;i<m.size();i++) m.at(i).setFileType(ZS);
      break;
    case uzs:
      for(i=0;i<m.size();i++) m.at(i).setFileType(UZS);
      break;
    case dunno:
    default:
      for(i=0;i<m.size();i++) m.at(i).setFileType(Unspecified);
      break;
  };
    
  switch(ff){
    case ms1:
    case ms2:
    case  zs:
    case uzs:
      setCompression(false);
      writeFile(c,true,m);
      break;
    case mzXML:
    case mzData:
      cout << "Cannot write mzXML or mzData formats. Nothing written." << endl;
      break;
    case bms1:
      setCompression(false);
      iFType=1;
      writeFile(c,false,m);
      break;
    case bms2:
      setCompression(false);
      iFType=3;
      writeFile(c,false,m);
      break;
    case cms1:
      setCompression(true);
      iFType=2;
      writeFile(c,false,m);
      break;
    case cms2:
      setCompression(true);
      iFType=4;
      writeFile(c,false,m);
      break;
    default:
      cout << "Unknown file format. Nothing written." << endl;
      break;
  };

};

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
	};
    
	fclose(fileOut);


};

void MSReader::appendFile(char* c, Spectrum& s){
  MSFileFormat ff;
	FILE* fileOut;

	if(c == NULL) return;
  ff=checkFileFormat(c);

  switch(ff){
    case ms1:
    case bms1:
    case cms1:
      s.setFileType(MS1);
      break;
    case ms2:
    case bms2:
    case cms2:
      s.setFileType(MS2);
      break;
    case zs:
      s.setFileType(ZS);
      break;
    case uzs:
      s.setFileType(UZS);
      break;
    case dunno:
    default:
      s.setFileType(Unspecified);
      break;
  };

  switch(ff){
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
    default:
      cout << "Cannot append file: unknown or unsupported file type." << endl;
      break;
  };

};

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
    case ms1:
    case ms2:
    case  zs:
    case uzs:
	    fileOut=fopen(c,"at");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,true,m.at(i));
        writeTextSpec(fileOut,m.at(i));
      };
      fclose(fileOut);
      break;
    case bms1:
    case bms2:
      fileOut=fopen(c,"ab");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,false,m.at(i));
        writeBinarySpec(fileOut,m.at(i));
      };
      fclose(fileOut);
      break;
    case cms1:
    case cms2:
      fileOut=fopen(c,"ab");
      for(i=0;i<m.size();i++){
        writeSpecHeader(fileOut,false,m.at(i));
        writeCompressSpec(fileOut,m.at(i));
      };
      fclose(fileOut);
      break;
    default:
      cout << "Cannot append file: unknown or unsupported file type." << endl;
      break;
  };

};

void MSReader::setPrecision(int i, int j){
	iIntensityPrecision=i;
	iMZPrecision=j;
};

void MSReader::setPrecisionInt(int i){
	iIntensityPrecision=i;
};

void MSReader::setPrecisionMZ(int i){
	iMZPrecision=i;
};

bool MSReader::readFile(char* c, Spectrum& s, int scNum){

  if(c!=NULL) lastFileFormat = checkFileFormat(c);   
  return readFile(c,lastFileFormat,s,scNum);

};

bool MSReader::readFile(char* c, MSFileFormat f, Spectrum& s, int scNum){

	//Redirect functions to appropriate places, if possible.
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
    case dunno:
		default:
			return false;
			break;
	};

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
      cout << "Error reading input file " << c << endl;
      return false;
		};
		rampFileOpen=true;
		//read the index
		indexOffset = getIndexOffset(rampFileIn);
		pScanIndex = readIndex(rampFileIn,indexOffset,&rampLastScan);
		rampIndex=0;

	} else {
		//if no new file requested, check to see if one is open already
		if (rampFileIn == NULL) return false;
	};


	//clear any spectrum data
	s.clear();

	//read scan header
	if(scNum!=0) {
		rampIndex=0;
		for(i=1;i<rampLastScan;i++){
			readHeader(rampFileIn, pScanIndex[i], &scanHeader);
			if(scanHeader.acquisitionNum==scNum) {
				rampIndex=i;
				break;
			};
		};
		if(rampIndex==0) return false;

		readHeader(rampFileIn, pScanIndex[rampIndex], &scanHeader);
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
		s.setScanNumber(scanHeader.acquisitionNum);
    s.setScanNumber(scanHeader.acquisitionNum,true);
		s.setRTime((float)scanHeader.retentionTime);
    if(strlen(scanHeader.activationMethod)>1){
      if(strcmp(scanHeader.activationMethod,"CID")==0) s.setActivationMethod(CID);
      if(strcmp(scanHeader.activationMethod,"ECD")==0) s.setActivationMethod(ECD);
      if(strcmp(scanHeader.activationMethod,"ETD")==0) s.setActivationMethod(ETD);
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

	} else {
		
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

			//if we got here, we passed the filter.
			break;
		};

		s.setScanNumber(scanHeader.acquisitionNum);
    s.setScanNumber(scanHeader.acquisitionNum,true);
		s.setRTime((float)scanHeader.retentionTime);
    if(strlen(scanHeader.activationMethod)>1){
      if(strcmp(scanHeader.activationMethod,"CID")==0) s.setActivationMethod(CID);
      if(strcmp(scanHeader.activationMethod,"ECD")==0) s.setActivationMethod(ECD);
      if(strcmp(scanHeader.activationMethod,"ETD")==0) s.setActivationMethod(ETD);
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

	};

	free(pPeaks);
	return true;

};

void MSReader::setFilter(MSSpectrumType m){
	filter=m;
};

void MSReader::setCompression(bool b){
	compressMe=b;
};

void MSReader::writeCompressSpec(FILE* fileOut, Spectrum& s){

	int j;

	//file compression
	int err;
	uLong len;
	Byte *comprM, *comprI;
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
	comprM = (Byte*)calloc((uInt)comprLenM, 1);
	err = compress(comprM, &comprLenM, (const Bytef*)pD, len);

	//compress intensity
	len = (uLong)s.size()*sizeof(float);
	sizeI = len;
	comprLenI = compressBound(len);
	comprI = (Byte*)calloc((uInt)comprLenI, 1);
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
	Byte *compr;
	double *mz;
	float *intensity;

	fread(&i,4,1,fileIn);
	mzLen = (uLong)i;
	fread(&i,4,1,fileIn);
	intensityLen = (uLong)i;
		
	compr = new Byte[mzLen];
	mz = new double[ms.numDataPoints];
	uncomprLen=ms.numDataPoints*sizeof(double);
	fread(compr,mzLen,1,fileIn);
	uncompress((Bytef*)mz, &uncomprLen, compr, mzLen);
	delete [] compr;
			
	compr = new Byte[intensityLen];
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

	int j,k;
	char t[64];

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
			};
		} else {
			fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
		};
	};
 
};

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
		mft=s.getFileType();
		if(mft==MS2 || mft==MS3 || mft==SRM){
			fprintf(fileOut,"S\t%d\t%d\t%.*f\n",s.getScanNumber(),s.getScanNumber(true),2,s.getMZ());
		} else {
			fprintf(fileOut,"S\t%d\t%d\n",s.getScanNumber(),s.getScanNumber(true));
		};
		if(s.getRTime()>0) fprintf(fileOut,"I\tRTime\t%.*f\n",4,s.getRTime());
		for(j=0;j<s.sizeZ();j++){
			fprintf(fileOut,"Z\t%d\t%.*f\n",s.atZ(j).z,2,s.atZ(j).mz);
		};
	} else {
    i=s.getScanNumber();
    fwrite(&i,4,1,fileOut);
    i=s.getScanNumber(true);
    fwrite(&i,4,1,fileOut);
    d=s.getMZ();
    fwrite(&d,8,1,fileOut);
    f=s.getRTime();
    fwrite(&f,4,1,fileOut);
    i=s.sizeZ();
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
		};
	};

};

void MSReader::readSpecHeader(FILE *fileIn, MSScanInfo &ms){

  fread(&ms.scanNumber[0],4,1,fileIn);
  if(feof(fileIn)) return;
  fread(&ms.scanNumber[1],4,1,fileIn);
  fread(&ms.mz,8,1,fileIn);
  fread(&ms.rTime,4,1,fileIn);
  fread(&ms.numZStates,4,1,fileIn);
  fread(&ms.numDataPoints,4,1,fileIn);

};

MSFileFormat MSReader::checkFileFormat(char *fn){

  int i,j;

  //check extension first - we must trust MS1 & MS2 & ZS & UZS
  i=strlen(fn);
  if(strcmp(fn+(i-4),".ms1")==0 || strcmp(fn+(i-4),".MS1")==0 ) return ms1;
  if(strcmp(fn+(i-4),".ms2")==0 || strcmp(fn+(i-4),".MS2")==0 ) return ms2;
  if(strcmp(fn+(i-3),".zs")==0 || strcmp(fn+(i-3),".ZS")==0 ) return zs;
  if(strcmp(fn+(i-4),".uzs")==0 || strcmp(fn+(i-4),".UZS")==0 ) return uzs;

  //For now, trust mzXML & mzData also
  if(strcmp(fn+(i-6),".mzXML")==0 || strcmp(fn+(i-6),".mzxml")==0 || strcmp(fn+(i-6),".MZXML")==0 ) return mzXML;
  if(strcmp(fn+(i-7),".mzData")==0 || strcmp(fn+(i-7),".mzdata")==0 || strcmp(fn+(i-7),".MZDATA")==0 ) return mzData;

  //We can check headers for other formats
  FILE *f;
  f=fopen(fn,"rb");
  if(f==NULL) {
    cout << "Error! Cannot open file to check file format." << endl;
    return dunno;
  };
  fread(&i,4,1,f);
  fread(&j,4,1,f);
  fclose(f);

  //cout << fn << "  FType = " << i << endl;
  switch(i){
    case 1: return bms1;
    case 2: return cms1;
    case 3: return bms2;
    case 4: return cms2;
    default: return dunno;
  };

  return dunno;

};
