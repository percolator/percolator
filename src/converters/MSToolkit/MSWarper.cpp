
/* 
use code from msconvertfile to produce output files 
1. args: take an msmat file initially, and then any number of ms1, ms2, etc.. files
   -will have to look up retention times for ms2 files from ms1 files, so argument structures will be:
   <msmat> <ms1> <ms2> - output modified versions of ms1, ms2, or have an in-place flag
*/




#include <iostream>
#include <string>
#include <iomanip>

#include "crawutils.H"
#include "msmat.H"
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"



using namespace std;



#ifdef _MSC_VER 
#ifndef MSCONVERT_MAIN
int mswarper_main(int argc, char * argv[] ) {
#else
int main(int argc, char *argv[]){
#endif
#else
int main(int argc, char *argv[]){
#endif






  MSFileFormat ff;

  MSReader w;
  MSObject o;
  Spectrum s;


  //How to use this program
  if(argc==1){
    printf("DESCRIPTION: Converts MS and MS/MS files from one MSToolkit format to another.\n\n");
    printf("USAGE: MSConvertFile <original file> <new file> <format>\n\n");
    printf("Valid formats are (case sensitive):\n\tt = text\n\tb = binary\n\tc = compressed\n");
    exit(1);
  }

  //Check to see if we have a valid input file
  //open the msmat file, without reading the actual intensity data into memory,
  //we just want the ort_wrt_map field

  crawutils::file_info mfile(argv[1]);
  msmat * m;
  msmat_open(&m, mfile, LM_SPARSE_DATA);
  if ( m->ort_wrt_map.size() == 0 ) {
    printf("msmat file %s needs to have a retention time map from an alignment\n",argv[1]);
    exit(1);
  }

  std::vector<const char *> ms1_file_names;

  for ( int i = 2 ; i < argc ; i++ ) {
    MSReader r;
    ff = r.checkFileFormat(argv[i]);
    if( ff==bms1 || ff==cms1 || ff==ms1 ) {
      printf("found MS1 file!\n");
      ms1_file_names.push_back(argv[i]);
    }
    else {
      printf("file %s must be MS1 file type\n",argv[i]);
      exit(0);
    }
  }
  //Set up the output file

  for ( int f_idx = 0; f_idx < ms1_file_names.size() ; f_idx++ ) {
    
    //since I don't know if this is reusable from file to file, keep it as a separate variable for each iteration
    MSReader r;

    ff = r.checkFileFormat((char*)ms1_file_names[f_idx]);
    r.readFile((char*)ms1_file_names[f_idx],s);

    string ms1_tmp_name(ms1_file_names[f_idx]);
    int last_dot_pos = ms1_tmp_name.rfind('.');
    ms1_tmp_name.insert(last_dot_pos,string(".warped"));
    char * ms1_name = (char*)ms1_tmp_name.c_str();
    
    /* find dot going left from rhs, then substitute in .warped before */

    if(s.getScanNumber()==0) {
      printf("Original file has no spectra!\n");
      exit(-2);
    }
    o.setHeader(r.getHeader());
    w.writeFile(ms1_name,ff,o);
    while(s.getScanNumber()!=0){
      w.appendFile(ms1_name,s);
      r.readFile(NULL,s);
    }
  }

  return 0;
}


