#include <iostream>
#include <iomanip>
#include <vector>
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace std;

const char * usage = "\n\
MSTAllRT : extracts Scan#,RT information from an ms2 file.\n\
\n\
Usage: MSTAllRT MS2_filename\n";

int main(int argc, char *argv[]){

  if ( argc != 2 ) {
    std::cerr << usage << std::endl;
    exit(1);
  }
	//Here are all the variable we are going to need
	MSReader r;
	Spectrum s;
	MSObject o;
	MSFileFormat ff;
	
	vector<float> rts(0);
	vector<int> scan_numbers(0);
	vector<double> mzs(0);
	ff = r.checkFileFormat(argv[1]);
	r.getHeader();
  r.readFile(argv[1],s);
  if(s.getScanNumber()==0) {
    printf("Original file has no spectra!\n");
    exit(-2);
  }
  rts.push_back(s.getRTime());
  scan_numbers.push_back(s.getScanNumber());
  mzs.push_back(s.getMZ());
    
  while(true) {
    r.readFile(NULL,s);
    if ( s.getScanNumber() == 0 ) {
        break;
     }
    rts.push_back(s.getRTime());
    scan_numbers.push_back(s.getScanNumber());
    mzs.push_back(s.getMZ());
  }
  for ( int i = 0 ; i < rts.size() ; i++ ) {
    std::cout << scan_numbers[i] <<"\t" << rts[i] << "\t" << mzs[i] << std::endl;
  }
  
  return 0;

};

  

