#include <string>
#include <iostream>
#include <map>
#include <fstream>

#include "DataSet.h"
#include "TmpDir.h"

class ValidateTabFile {       
  public:     
    std::string concatenateMultiplePINs(std::vector<std::basic_string<char>> fileNames);        
};