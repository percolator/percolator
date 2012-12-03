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
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"
#include "Caller.h"
using namespace std;


int main(int argc, char** argv) {
  
  Caller* pCaller = new Caller();
  int retVal = -1;
  
  try
  {
    if (pCaller->parseOptions(argc, argv)) {
      retVal = pCaller->run();
    }
  }
  catch (const std::exception& e) 
  {
    std::cerr << e.what() << endl;
    retVal = -1;
  }
  catch(...)
  {
    std::cerr << "Unknown exception, contact the developer.." << std::endl;
    retVal = -1;
  }
  
  delete pCaller;
  Globals::clean();
  return retVal;
}
