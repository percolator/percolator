#include <iostream>
#include <iomanip>
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace std;

//#ifdef _MSC_VER 
//#ifndef MSSINGLESCAN_MAIN
//int mssinglescan_main(int argc, char * argv[] ) {
//#else
//int main(int argc, char *argv[]){
//#endif
//#else
int main(int argc, char *argv[]){
//#endif


	//Here are all the variable we are going to need
	MSReader r;
	Spectrum s;
  int j;

  if(argc==1){
    printf("DESCRIPTION: Reads an MS/MS spectrum from any MSToolkit supported file type and outputs to screen in MS2 format.\n\n");
    printf("USAGE: MSSingleScan [scan number] [file]\n");
    exit(0);
  };

	r.readFile(argv[2],s,atoi(argv[1]));
  if(s.getScanNumber()==0) exit(-1);
  printf("S\t%d\t%d\t%.*f\n",s.getScanNumber(),s.getScanNumber(),2,s.getMZ());
	if(s.getRTime()>0) printf("I\tRTime\t%.*f\n",4,s.getRTime());
	for(j=0;j<s.sizeZ();j++){
		printf("Z\t%d\t%.*f\n",s.atZ(j).z,2,s.atZ(j).mz);
	};

	for(j=0;j<s.size();j++){
		printf("%.4f %.4f\n",s.at(j).mz,s.at(j).intensity);
	};

  return 0;

};
  

