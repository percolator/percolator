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
#include <cstdlib>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include "Caller.h"
#include "Globals.h"
using namespace std;


int main(int argc, char** argv) {
#ifdef __GLIBC__
  // glibc serves a block above its mmap threshold from a mapping of its own
  // and returns it to the system when it is freed. The threshold starts at
  // 128 KiB but rises to the size of every such block that is freed, up to
  // 32 MiB, and smaller blocks then stay in the heap after they are freed:
  // the training data of the first stages remained resident during the
  // later ones. A fixed threshold switches the adjustment off. Set here and
  // not in the library, whose callers choose their own allocator settings.
  mallopt(M_MMAP_THRESHOLD, 1024 * 1024);
#endif
  int retVal = EXIT_FAILURE; 
  {
    Caller caller;
    try
    {
      if (caller.parseOptions(argc, argv)) {
        if (caller.run()) retVal = EXIT_SUCCESS;
      }
    }
    catch (const std::exception& e) 
    {
      std::cerr << "Exception caught: " << e.what() << endl;
      retVal = EXIT_FAILURE;
    }
    catch(...)
    {
      std::cerr << "Unknown exception, contact the developer.." << std::endl;
      retVal = EXIT_FAILURE;
    }
  }
  Globals::clean();
  return retVal;
}
