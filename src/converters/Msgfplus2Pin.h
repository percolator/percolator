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

#ifndef MSGFPLUS2PIN_H_
#define MSGFPLUS2PIN_H_

#include "MsgfplusReader.h"
#include "Interface.h"

using namespace std;

class Msgfplus2pin : public Interface {

 public:
   
	Msgfplus2pin();
	virtual ~Msgfplus2pin();
	static std::string Usage()
	{
	  std::stringstream endnote;
	  endnote << "Usage:" << endl;
	  endnote << "   msgf2pin [options] target.mzid decoy.mzid" << endl << endl;
	  endnote << "target.mzid and decoy.mzid are MzIdentML-files of MS-GF+ from" << endl;
	  endnote << "separate target and decoy searches. Internal MS-GF+ target/decoy" << endl;
	  endnote << "analysis should be turned off, and the addFeatures options turned on." << endl;
	  endnote << "Multiple MzIdentML-files can be merged by replacing target and decoy" << endl;
	  endnote << "filepaths with meta files. Meta files are text files containing the" << endl;
	  endnote << "the paths of mzid-files, one path per line. For successful results," << endl;
	  endnote << "the different runs should be generated under similar conditions." << endl;
	  return endnote.str();
	}
	int run();
	
 private:
  
	MsfgplusReader *reader;
};

int main(int argc, char **argv);

#endif /* MSGFPLUS2PIN_H_ */

