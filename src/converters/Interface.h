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

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include <string>
#include "Option.h"
#include "config.h"
#include "parseoptions.h"
#include "Enzyme.h"
#include "Globals.h"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
using namespace std;

class Interface {

 public:
   
	Interface();
	virtual ~Interface();
	virtual std::string greeter();
	virtual std::string extendedGreeter();
	virtual bool parseOpt(int argc, char **argv,const std::string &usage);
	
 protected:
  
	ParseOptions parseOptions;
	std::string targetFN;
	std::string decoyFN;
	std::string xmlOutputFN;
	std::string call;
	std::string spectrumFile;
};

#endif /* INTERFACE_H_ */
