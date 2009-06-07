#include <iostream>
#include <iomanip>
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace std;

//MH: What is msconvert_main? I cannot find any docs on it...
//#define MSCONVERT_MAIN

//#ifdef _MSC_VER 
//#ifndef MSCONVERT_MAIN
//int msconvert_main(int argc, char * argv[] ) {
//#else
//int main(int argc, char *argv[]){
//#endif
//#else
int main(int argc, char *argv[]){
//#endif

	//Here are all the variable we are going to need
  MSFileFormat ff;
	MSReader r;
  MSReader w;
  MSObject o;
	Spectrum s;
  int j;

  //How to use this program
  if(argc==1){
    printf("DESCRIPTION: Converts MS and MS/MS files from one MSToolkit format to another.\n\n");
    printf("USAGE: MSConvertFile <original file> <new file> <format>\n\n");
    printf("Valid formats are (case sensitive):\n\tt = text\n\tb = binary\n\tc = compressed\n");
    exit(0);
  }

  //Check to see if we have a valid input file
  ff = r.checkFileFormat(argv[1]);
  if(ff==dunno) {
    printf("Unknown file format!\n");
    exit(-1);
  }

  //Check to see if we have a valid input file
//  if(ff==mzXML || ff==mzData || ff==zs || ff==uzs) {
  if(ff==mzData || ff==zs || ff==uzs) {
    printf("File format not supported. Please use .ms1, .ms2, .bms1, .bms2, .cms1, .cms2\n");
    exit(-4);
  }

  //Check if format is valid
  j=0;
  if(argv[3][0]=='b') j=1;
  if(argv[3][0]=='c') j=2;
  if(argv[3][0]=='t') j=3;
  if(j==0){
    printf("Invalid format requested.\n");
    printf("Valid formats are (case sensitive):\n\tt = text\n\tb = binary\n\tc = compressed\n");
    exit(-3);
  }

  //Set up the output file
  std::vector<MSSpectrumType> filts;
  filts.push_back(MS1);
  r.setFilter(filts);
	r.readFile(argv[1],s);
  if(s.getScanNumber()==0) {
    printf("Original file has no spectra!\n");
    exit(-2);
  }
  o.setHeader(r.getHeader());
  switch(ff){
    case ms1:
    case bms1:
    case cms1:
      if(j==1) w.writeFile(argv[2],bms1,o);
      if(j==2) w.writeFile(argv[2],cms1,o);
      if(j==3) w.writeFile(argv[2],ms1,o);
      break;
    case ms2:
    case bms2:
    case cms2:
      if(j==1) w.writeFile(argv[2],bms2,o);
      if(j==2) w.writeFile(argv[2],cms2,o);
      if(j==3) w.writeFile(argv[2],ms2,o);
      break;
    default:
      break;
  }

  while(s.getScanNumber()!=0){
    w.appendFile(argv[2],s);
    r.readFile(NULL,s);
  }

  return 0;

};
  

