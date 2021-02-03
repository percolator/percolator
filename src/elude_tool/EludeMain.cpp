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
#include "EludeCaller.h"

int main(int argc, char** argv) {
  
  EludeCaller* eCaller = new EludeCaller();
  int retVal = -1;
  
  try
  {
    if (eCaller->ParseOptions(argc, argv)) {
      retVal = eCaller->Run();
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

  delete eCaller;
  Globals::clean();
  return retVal;
}

