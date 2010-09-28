#include <iostream>
#include <iomanip>
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace std;


//#ifdef _MSC_VER 
//#ifndef MSTDEMO_MAIN
//int msconvert_main(int argc, char * argv[] ) {
//#else
//int main(int argc, char *argv[]){
//#endif
//#else
int main(int argc, char *argv[]){
//#endif

	//Here are all the variable we are going to need
	MSReader r;
	Spectrum s;
	MSObject o;

  /*
  r.readFile(argv[1],s);
  o.setHeader(r.getHeader());
  MSReader r2;
  r2.writeFile("test.bms1",bms1,o);
  
  while(s.getScanNumber()!=0){
    r2.appendFile("test.bms1",s);
    r.readFile(NULL,s);
  };

  r.readFile(argv[1],s);
  o.setHeader(r.getHeader());
  r2.writeFile("test.cms1",cms1,o);
  
  while(s.getScanNumber()!=0){
    r2.appendFile("test.cms1",s);
    r.readFile(NULL,s);
  };
  */

  
  int i;
  r.readFile("test.cms1",s,770);
  cout << s.getScanNumber() << "\t" << s.getRTime() << endl;
  for(i=0;i<s.size();i++) cout << s[i].mz << " " << s[i].intensity << endl;
  cout << "\n\n" << endl;

  /*
  r.readFile(argv[1],s,770);
  cout << s.getScanNumber() << "\t" << s.getRTime() << endl;
  for(i=0;i<s.size();i++) cout << s[i].mz << " " << s[i].intensity << endl;
  cout << "\n\n" << endl;

  r.readFile("test.bms1",s,770);
  cout << s.getScanNumber() << "\t" << s.getRTime() << endl;
  for(i=0;i<s.size();i++) cout << s[i].mz << " " << s[i].intensity << endl;
  cout << "\n\n" << endl;

  r.readFile("test.cms1",s,770);
  cout << s.getScanNumber() << "\t" << s.getRTime() << endl;
  for(i=0;i<s.size();i++) cout << s[i].mz << " " << s[i].intensity << endl;
  cout << "\n\n" << endl;
  */

	exit(1);

	//read the first scan
	cout << "Read first scan of " << argv[1] << endl;
	if(!r.readFile(argv[1],true,s)) cout << "Error reading scan" << endl;

	//output the number of data points
	cout << "Scan #" << s.getScanNumber() << " has " << s.size() << " data points." << endl;

	//add header to our object
	cout << "Copy file header to our object...";
	o.setHeader(r.getHeader());
	cout << "Done!" << endl;

	//add spectrum to our MSObject
	cout << "Adding spectrum to our object...";
	o.add(s);
  cout << "Done!" << endl;

	//Validate spectrum
	cout << "Spectrum in MSObject is scan #" << o.at(0).getScanNumber() 
			 << " and has " << o.at(0).size() << " data points." << endl;

	//Add next 10 spectra
	cout << "Load next 10 spectra:" << endl;
	for(int i=0;i<10;i++){
		if(!r.readFile(NULL,true,s)) cout << "Error reading scan" << endl;
		cout << "  loaded scan #" << s.getScanNumber() << "...";
		o.add(s);
		cout << "added to MSObject!" << endl;
	};

	//Check 5th scan of MSObject
	cout << "The 5th spectrum in MSObject is scan #" << o.at(4).getScanNumber() 
			 << " and has " << o.at(4).size() << " data points." << endl;

	//Set intensity precision to 4 decimal places
	cout << "Setting intensity precision to 4 decimal places...";
	r.setPrecisionInt(4);
	cout << "Done!" << endl;

	//Write file as text
	cout << "Writing file to out.txt as text...";
	r.writeFile("out.txt",true,o);
	cout << "Done!" << endl;

	//Write file as binary
	cout << "Writing file to out.bin as binary...";
	r.writeFile("out.bin",false,o);
	cout << "Done!" << endl;

	//Read in one more spectrum, store to temporary variable
	Spectrum tempSpec;
	cout << "Reading 12th spectrum...";
	if(!r.readFile(NULL,true,tempSpec)) cout << "Error reading scan" << endl;
	cout << "Done! Stored to memory." << endl;

	//Read in binary file and check scan numbers.
	cout << "Clearing object and reading in entire binary file...";
	o.clear();
	if(!r.readFile("out.bin",false,s)) cout << "Error reading scan" << endl;
	o.add(s);
	while(r.readFile(NULL,false,s)){
		o.add(s);
	};
	cout << "Done!" << endl;

	//Check the 5th scan
	cout << "The 5th spectrum in the binary file is scan #" << o.at(4).getScanNumber() 
			 << " and has " << o.at(4).size() << " data points." << endl;

	//Append our 12th spectrum to the end of the binary file
	cout << "Appending 12th spectrum to the end of the binary file...";
	r.appendFile("out.bin",false,tempSpec);
	cout << "Done!" << endl;

	cout << "This ends the tutorial. BLARG!" << endl;

  return 0;

};



