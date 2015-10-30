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

#include "Option.h"
#include "Globals.h"
#include "FisherCaller.h"

int main(int argc, char** argv) {
  std::map<std::string, std::string> fragment_map, duplicate_map;
  
  FisherCaller* fCaller = new FisherCaller(TRYPSIN, FULL_DIGEST, 10, 30, 0);
  int retVal = EXIT_FAILURE;
  
  try
  {
    if (fCaller->parseOptions(argc, argv)) {
      retVal = fCaller->getProteinFragmentsAndDuplicates(fragment_map, duplicate_map);
    }
  }
  catch (const std::exception& e) 
  {
    std::cerr << e.what() << endl;
    retVal = EXIT_FAILURE;
  }
  catch(...)
  {
    std::cerr << "Unknown exception, contact the developer.." << std::endl;
    retVal = EXIT_FAILURE;
  } 

  delete fCaller;
  return retVal;
}

