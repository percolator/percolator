/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef SQT2PIN_H_
#define SQT2PIN_H_

#include "SqtReader.h"
#include "Interface.h"

using namespace std;

class Sqt2Pin : public Interface {

 public:
   
	Sqt2Pin();
	virtual ~Sqt2Pin();
	static std::string Usage()
	{
	  ostringstream endnote;
	  endnote << "Usage:" << endl;
	  endnote << "   sqt2pin [options] -o output.xml target.sqt decoy.sqt " << endl << endl;
	  endnote << "Where output.xml is where the output will be written (ensure to have read and " << endl;
	  endnote << "write access on the file).target.sqt is the target sqt-file, and decoy.sqt is" << endl;
	  endnote << "the decoy sqt-file. Small data sets may be merged by replace the sqt-files with" << endl;
	  endnote << "meta files. Meta files are text files containing the paths of sqt-files, one" << endl;
	  endnote << "path per line. For successful result, the different runs should be generated" << endl;
	  endnote << "under similar condition." << endl;
	  return endnote.str();
	}
	int run();
	
 private:
  
	SqtReader *reader;
};

int main(int argc, char **argv);

#endif /* SQT2PIN_H_ */
