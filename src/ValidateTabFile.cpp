#include "ValidateTabFile.h"

/**
 * Concatenates multiple tab-files into one file in outputDir, using the
 * header of the file with the most columns. Columns missing in a file (e.g.
 * charge states) are filled with zeros.
 * @return name of the concatenated file
 */
std::string ValidateTabFile::concatenateMultiplePINs(
    const std::vector<std::string>& fileNames, const std::string& outputDir) {
  // Validate multiple tab-files, and then concatenate them.
  std::string inputFN;
  std::ifstream pinFileStream;

  /* Complete column position to header name */
  std::map<int,std::string> columnMap = {};
  /* Temporary mapper while searching for complete position-to-header mapper */
  std::map<int,std::string> tmpMap;

  /* Keep track on col index */
  int column_index;

  /* Loop through all files to search for column headers */
  for(const std::string &text : fileNames) {
      
    /* Get header */
    pinFileStream.open(text.c_str(), std::ios::in);
    std::string headerRow;

    getline(pinFileStream, headerRow);

    TabReader reader(headerRow);

    column_index = 0;
    while (!reader.error()) {
      std::string optionalHeader = reader.readString();

      tmpMap[column_index] = optionalHeader;
      column_index++;
    }

    /* Get header from file with the most columns */
    if (columnMap.size() < tmpMap.size()) {
      columnMap = tmpMap;
    }
    tmpMap.clear();
    pinFileStream.close();
  }

  /* Create a temporary tab file to concatenate tab files in to */ 
  inputFN = outputDir + "/converters-tmp.tcb";
  ofstream outFile;
  outFile.open(inputFN);
  if (!outFile.is_open()) {
    throw MyException("ERROR: Could not create the temporary file " + inputFN +
                      " for concatenating the input files.\n");
  }

  /* Only want to print header once */
  bool firstFile = true;
  /* If go to next column when going though missing charge states */
  bool nextColumn;

  /* Start appending pin files to one pin file. */
  for(const string &text : fileNames) {
    if (VERB > 0) {
      std::cerr << "Reading file: " << text.c_str() << std::endl;
    }
    /* Keep track of missing charge columns for a given file */
    std::vector<int> missingCols;
        
    /* Read pin file */
    pinFileStream.open(text.c_str(), std::ios::in);
    std::string headerRow;
    getline(pinFileStream, headerRow);
    TabReader reader(headerRow);

    /* Keep track on column index */
    column_index = 0;
    /* Go through all columns to print header and search for missing charge state */
    while (!reader.error()) {
      std::string optionalHeader = reader.readString();
          
      /* Check if charge state is missing */
      if (columnMap[column_index] != optionalHeader) {
        /* Charge state is missing at  'column_index' */
        missingCols.push_back(column_index);
        nextColumn=true;
        if (firstFile) {
          outFile << columnMap[column_index] << "\t";
        }
        /* Check if there are any adjacent missing charge states */
        while(nextColumn) {
          column_index++;
          if (column_index >= static_cast<int>(columnMap.size())) {
            throw MyException("ERROR: The column " + optionalHeader + " of " +
                text + " is not in the header with the most columns, " +
                "the input files cannot be concatenated.\n");
          }
          if (columnMap[column_index] == optionalHeader) {
            /* No missing charge state, continue */
            nextColumn=false;
          } else {
            /* Missing charge state */
            missingCols.push_back(column_index);
            if (firstFile) {
              outFile << columnMap[column_index] << "\t";
            }             
          }
        }
      } 
      if (firstFile) {
        outFile << columnMap[column_index] << "\t";
      }
      column_index++;
    }

    /* Go through rest of the rows in pin-file */
    std::string nextRow;
    getline(pinFileStream, nextRow);

    if (firstFile) {
      outFile <<  "\n";
    } else {
      std::string defaultDirectionString = "defaultdirection";
      if (nextRow.size() >= defaultDirectionString.size()) {
        std::string psmid = nextRow.substr(0, defaultDirectionString.size());
        std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
        if (psmid == defaultDirectionString) {
          getline(pinFileStream, nextRow); // skip row with default direction
        }
      }
    }

    do {
      if (missingCols.empty()) {
        // same as splitting at the tabs and writing each field followed by
        // a tab, as done below
        outFile << nextRow << "\t\n";
        continue;
      }
      TabReader readerRow(nextRow);

      /* Keep track on column for missing charge states */
      int col = 0;
      while (!readerRow.error()) {           
        std::string value = readerRow.readString();
        /* Check for missing charge state */
        if (std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
          nextColumn=true;
          outFile << 0 << "\t";
          col++;
          /* Check for adjecent missing charge states */
          while(nextColumn) {
            if(std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
              /* Missing adjecent charge state  */
              outFile << 0 << "\t";
              col++;
            } else {
              /* Charge states  */
              nextColumn = false;
            }
          }  
        }
        outFile << value << "\t";
        col++;
      }
      outFile << "\n";
    } while (getline(pinFileStream, nextRow));
    pinFileStream.close();
    firstFile = false;      
  }
  outFile.close();
  if (!outFile) {
    throw MyException("ERROR: Could not write the temporary file " + inputFN +
                      ".\n");
  }
  return inputFN;
}
