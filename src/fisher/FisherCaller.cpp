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
#include <cstdio>
#include <iostream>

#include "Database.h"
using namespace std;

int main(int argc, char** argv) {
  
  // now create a database
  Database* database_ = new Database(argv[1], false);
  
  if(!database_->parse()){
    std::cerr <<  "Failed to parse database, cannot create new index for " << argv[1] << std::endl;
  } else {
    std::cerr << "Read " << database_->getNumProteins() << " proteins." << std::endl;
    std::cerr << database_->getProteinAtIdx(5)->getSequencePointer() << std::endl;
  }
  
  return EXIT_SUCCESS;
}
