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

#ifndef MZIDENTML2PIN_H
#define MZIDENTML2PIN_H


#include "MzidentmlReader.h"
#include "Option.h"
#include "config.h"
#include <string>

using namespace std;

class Mzidentml2pin {
  

 public:
	Mzidentml2pin();
	
	virtual ~Mzidentml2pin();
	
	std::string greeter();
	
	std::string extendedGreeter();
	
	bool parseOpt(int argc, char **argv);
	
	int run();


 private:
  
	ParseOptions parseOptions;
	std::string targetFN;
	std::string decoyFN;
	std::string xmlOutputFN;
	std::string call;
	std::string spectrumFile;
	MzidentmlReader *reader;
};

int main(int argc, char **argv);

#endif /*MZIDENTML2PIN_H*/
