/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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

#ifndef TANDEM2PIN_H_
#define TANDEM2PIN_H_

#include "TandemReader.h"
#include "Interface.h"
#include <string>

using namespace std;

class Tandem2Pin : public Interface {

  public:
    
	Tandem2Pin();
	virtual ~Tandem2Pin();

	static std::string Usage()
	{
	  ostringstream endnote;
	  endnote << "Usage:" << endl;
	  endnote << "   tandem2pin [options] -o output.xml target_file decoy_file " << endl << endl;
	  endnote << "Where output.xml is where the output will be written (ensure to have read and " << endl;
	  endnote << "write access on the file).target_file is the target X!tandem-file, and decoy_file is" << endl;
	  endnote << "the decoy X!tandem-file. Small data sets may be merged by replace the X!tandem-files with" << endl;
	  endnote << "meta files. Meta files are text files containing the paths of X!tandem-files, one" << endl;
	  endnote << "path per line. For successful result, the different runs should be generated" << endl;
	  endnote << "under similar condition." << endl;
	  return endnote.str();
	}

	int run();

  private:

	TandemReader *reader;
};

int main(int argc, char **argv);

#endif /* TANDEM2PIN_H_ */
