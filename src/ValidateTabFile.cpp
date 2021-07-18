#include "ValidateTabFile.h"

std::string ValidateTabFile::concatenateMultiplePINs(std::vector<std::basic_string<char>> fileNames) {

  std::string inputFN_;
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
    std::string hederRow;

    getline(pinFileStream, hederRow);

    TabReader reader(hederRow);

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

  string tcf = "";
  char tcd;

  TmpDir tmpDir;
  tmpDir.createTempFile(&tcf, &tcd);
  inputFN_ = tcf;
  ofstream outFile;
  outFile.open(inputFN_);

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
    std::string hederRow;
    getline(pinFileStream, hederRow);
    TabReader reader(hederRow);
        
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

    if (firstFile) {
      outFile <<  "\n";
    }
    /* Go through rest of the rows in pin-file */
    std::string nextRow;
        
    while(getline(pinFileStream, nextRow)) {
      
      TabReader readerRow(nextRow);

      /* Keep track on column for missing charge states */
      int col = 0;
      while (!readerRow.error()) {           
        std::string value = readerRow.readString();
        /* Check for missing charge state */
        if(std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
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
    }
    pinFileStream.close();
    firstFile = false;      
  }
  outFile.close();
  return inputFN_;
}
