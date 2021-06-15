#include <boost/filesystem.hpp>
#include <string>
#include <iostream>

class TmpDir {       
  public:             
    void createTempFile(std::string* tcf, char* tcd);
};