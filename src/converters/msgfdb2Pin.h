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

#ifndef MSGFDB2PIN_H_
#define MSGFDB2PIN_H_

#include "msgfdbReader.h"
#include "Option.h"
#include "config.h"

using namespace std;

class msgfdb2Pin {
public:
	msgfdb2Pin();
	virtual ~msgfdb2Pin();
	std::string greeter();
	std::string extendedGreeter();
	bool parseOpt(int argc, char **argv);
	int run();

protected:
	ParseOptions parseOptions;
	std::string targetFN;
	std::string decoyFN;
	std::string xmlOutputFN;
	std::string call;
	std::string spectrumFile;
	msgfdbReader *reader;
};

int main(int argc, char **argv);

#endif /* MSGFDB2PIN_H_ */
