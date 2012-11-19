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

#ifndef SEQUEST2PIN_H_
#define SEQUEST2PIN_H_

#include "SequestReader.h"
#include "Interface.h"

using namespace std;

class Sequest2Pin : public Interface {

 public:
   
	Sequest2Pin();
	virtual ~Sequest2Pin();
	static std::string Usage()
	{
	  ostringstream endnote;
	  endnote << "Usage:" << endl;
	  endnote << "   sequest2pin [options] target.mzid decoy.mzid" << endl << endl;
	  endnote << "Target.mzid [Sequest] is the target MzIdentML-file, and decoy.mzid [Sequest] is" << endl;
	  endnote << "the decoy MzIdentML-file. Small data sets may be merged by replace the mzid-files with" << endl;
	  endnote << "meta files. Meta files are text files containing the paths of mzid-files, one" << endl;
	  endnote << "path per line. For successful result, the different runs should be generated" << endl;
	  endnote << "under similar conditions." << endl;
	  return endnote.str();
	}
	int run();
	
 private:
  
	SequestReader *reader;
};

int main(int argc, char **argv);

#endif /* SEQUEST2PIN_H_ */
