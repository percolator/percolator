#include <boost/filesystem.hpp>
#include <string>
#include <iostream>
#include <map>
#include "DataSet.h"
#include <fstream>

class TmpDir {       
  public:             
    void static createTempFile(std::string* tcf, char* tcd);
    std::string concatenateMultiplePINs(std::vector<std::basic_string<char>> fileNames);
};
