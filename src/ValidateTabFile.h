#ifndef VALIDATETABFILE_H_
#define VALIDATETABFILE_H_

#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <vector>

#include "DataSet.h"

// The lines of one input file in a concatenated input: line firstLine of the
// concatenated file is line firstFileLine of fileName, and so on.
struct InputSegment {
  unsigned int firstLine;
  unsigned int firstFileLine;
  std::string fileName;
};

class ValidateTabFile {
  public:
    // If segments is given, it receives one entry per input file with PSM
    // lines, in the order of the concatenated file.
    std::string concatenateMultiplePINs(const std::vector<std::string>& fileNames,
                                        const std::string& outputDir,
                                        std::vector<InputSegment>* segments = nullptr);
};

#endif /*VALIDATETABFILE_H_*/
