#include <string>
#include <iostream>
#include <map>
#include <fstream>

#include "DataSet.h"

class ValidateTabFile {       
  public:     
    std::string concatenateMultiplePINs(const std::vector<std::string>& fileNames,
                                        const std::string& outputDir);
};