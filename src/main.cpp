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
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"
#include "Caller.h"

int main(int argc, char** argv) {
  Caller* pCaller = new Caller(false);
  Caller* pCallerPeptide = new Caller(true);
  int retVal = -1;
  // reading command line arguments
  bool validArguments = pCaller->parseOptions(argc, argv);
  validArguments = validArguments && pCallerPeptide->parseOptions(argc, argv);
  if (validArguments) {
    // executing: psm
    retVal = pCaller->run();
    Globals::clean();
    // executing: unique peptides
    retVal += pCallerPeptide->run();
  }
  // outputing results to XML file
  if (pCallerPeptide->xmloutFN.size() > 0) {
    ofstream xmlStream(pCallerPeptide->xmloutFN.data(), ios::out);
    pCallerPeptide->writeXML(xmlStream, pCaller->fullset, pCallerPeptide->fullset);
    xmlStream.close();
  }
  delete pCaller;
  delete pCallerPeptide;
  Globals::clean();
  return retVal;
}
