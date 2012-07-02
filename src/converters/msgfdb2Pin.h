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

#ifndef PIN_VERSION_MAJOR
#define PIN_VERSION_MAJOR "@PIN_VERSION_MAJOR@"
#endif
#ifndef PIN_VERSION_MINOR
#define PIN_VERSION_MINOR "@PIN_VERSION_MINOR@"
#endif
#ifndef WRITABLE_DIR
#define WRITABLE_DIR "@WRITABLE_DIR@"
#endif
#ifndef TEMP_DIR
#define TEMP_DIR "@TEMP_DIR@"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include "msgfdbReader.h"
#include "Option.h"
#include "config.h"
#include "serializer.hxx"
#include "MSReader.h"
#include "Spectrum.h"
#include "MSToolkitTypes.h"
#include "DataSet.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include <cmath>

using namespace std;

class msgfdb2Pin {
public:
	msgfdb2Pin();
	virtual ~msgfdb2Pin();
	std::string greeter();
	std::string extendedGreeter();
	bool parseOpt(int argc, char **argv);
	int run();
	void readRetentionTime(std::string filename);
	void storeRetentionTime(FragSpectrumScanDatabase* database);

protected:
	ParseOptions parseOptions;
	std::string targetFN;
	std::string decoyFN;
	std::string xmlOutputFN;
	std::string call;
	std::string spectrumFile;
	map<int, vector<double> > scan2rt;
	msgfdbReader *reader;
};

int main(int argc, char **argv);

#endif /* MSGFDB2PIN_H_ */
