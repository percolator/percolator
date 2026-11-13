#include "ValidateTabFile.h"

#include <algorithm>

/**
 * Reads the headers of the files: the header of the file with the most
 * columns is the header of the concatenated input, and the columns that a
 * file lacks (e.g. charge states) are filled with zeros in its rows. Each
 * file is named on stderr, a file without PSMs gets a warning, and a file
 * with a column that the widest header lacks is an error, as when the files
 * were copied into one file before reading.
 */
ConcatenatedPins::ConcatenatedPins(const std::vector<std::string>& fileNames) {
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

  /* Only want to print header once */
  bool firstFile = true;
  /* If go to next column when going though missing charge states */
  bool nextColumn;

  /* The header, and what to do with the rows of each file. */
  for(const string &text : fileNames) {
    if (VERB > 0) {
      std::cerr << "Reading file: " << text.c_str() << std::endl;
    }
    File file;
    file.name = text;
    /* Keep track of missing charge columns for a given file */
    std::vector<int>& missingCols = file.missingCols;
        
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
          headerLine_ += columnMap[column_index] + "\t";
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
              headerLine_ += columnMap[column_index] + "\t";
            }             
          }
        }
      } 
      if (firstFile) {
        headerLine_ += columnMap[column_index] + "\t";
      }
      column_index++;
    }

    /* Look at the first rows of the pin-file */
    std::string nextRow;
    bool hasRow = static_cast<bool>(getline(pinFileStream, nextRow));
    unsigned int fileLine = 2;

    if (firstFile) {
      headerLine_ += "\n";
    } else {
      std::string defaultDirectionString = "defaultdirection";
      if (nextRow.size() >= defaultDirectionString.size()) {
        std::string psmid = nextRow.substr(0, defaultDirectionString.size());
        std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
        if (psmid == defaultDirectionString) {
          // skip row with default direction
          hasRow = static_cast<bool>(getline(pinFileStream, nextRow));
          ++fileLine;
        }
      }
    }
    pinFileStream.close();
    firstFile = false;
    // a file with only a header would otherwise add an empty PSM line
    if (!hasRow) {
      std::cerr << "Warning: the input file " << text
                << " contains no PSMs." << std::endl;
      continue;
    }
    file.firstRowLine = fileLine;
    files_.push_back(file);
  }
}

void ConcatenatedPins::rewind() {
  current_.close();
  current_.clear();
  nextFile_ = 0u;
  outLines_ = 0u;
  headerDone_ = false;
  fileOpen_ = false;
}

/*
 * Opens the next file with rows and moves to its first row. Returns false
 * after the last file.
 */
bool ConcatenatedPins::openNextFile() {
  current_.close();
  current_.clear();
  fileOpen_ = false;
  if (nextFile_ >= files_.size()) {
    return false;
  }
  const File& file = files_[nextFile_];
  current_.open(file.name.c_str(), std::ios::in);
  std::string skipped;
  for (unsigned int fileLine = 1; fileLine < file.firstRowLine; ++fileLine) {
    getline(current_, skipped);
  }
  // the lines of this file start here; known for all files after one pass
  if (segments_.size() == nextFile_) {
    segments_.push_back(
        InputSegment{outLines_ + 1, file.firstRowLine, file.name});
  }
  ++nextFile_;
  fileOpen_ = true;
  return true;
}

bool ConcatenatedPins::appendNextLine(std::string& out) {
  if (!headerDone_) {
    headerDone_ = true;
    out += headerLine_;
    ++outLines_;
    return true;
  }
  while (!fileOpen_ || !getline(current_, row_)) {
    if (!openNextFile()) {
      return false;
    }
  }
  const std::vector<int>& missingCols = files_[nextFile_ - 1].missingCols;
  ++outLines_;
  if (missingCols.empty()) {
    // same as splitting at the tabs and writing each field followed by
    // a tab, as done below
    out += row_;
    out += "\t\n";
    return true;
  }
  TabReader readerRow(row_);

  /* Keep track on column for missing charge states */
  int col = 0;
  bool nextColumn;
  while (!readerRow.error()) {           
    std::string_view value = readerRow.readView();
    /* Check for missing charge state */
    if (std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
      nextColumn=true;
      out += "0\t";
      col++;
      /* Check for adjecent missing charge states */
      while(nextColumn) {
        if(std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
          /* Missing adjecent charge state  */
          out += "0\t";
          col++;
        } else {
          /* Charge states  */
          nextColumn = false;
        }
      }  
    }
    out.append(value.data(), value.size());
    out += '\t';
    col++;
  }
  out += '\n';
  return true;
}

ConcatenatedPinStream::Buffer::int_type
ConcatenatedPinStream::Buffer::underflow() {
  if (gptr() < egptr()) {
    return traits_type::to_int_type(*gptr());
  }
  text_.clear();
  while (text_.size() < kBlockSize && pins_.appendNextLine(text_)) {
  }
  if (text_.empty()) {
    return traits_type::eof();
  }
  setg(&text_[0], &text_[0], &text_[0] + text_.size());
  return traits_type::to_int_type(*gptr());
}

ConcatenatedPinStream::Buffer::pos_type
ConcatenatedPinStream::Buffer::seekoff(off_type off,
                                       std::ios_base::seekdir dir,
                                       std::ios_base::openmode) {
  // only back to the start, for the second pass of -N
  if (off != 0 || dir != std::ios_base::beg) {
    return pos_type(off_type(-1));
  }
  pins_.rewind();
  setg(NULL, NULL, NULL);
  return pos_type(off_type(0));
}

ConcatenatedPinStream::Buffer::pos_type
ConcatenatedPinStream::Buffer::seekpos(pos_type pos,
                                       std::ios_base::openmode which) {
  return seekoff(off_type(pos), std::ios_base::beg, which);
}

/**
 * Concatenates multiple tab-files into one file in outputDir, using the
 * header of the file with the most columns. Columns missing in a file (e.g.
 * charge states) are filled with zeros.
 * @return name of the concatenated file
 */
std::string ValidateTabFile::concatenateMultiplePINs(
    const std::vector<std::string>& fileNames, const std::string& outputDir,
    std::vector<InputSegment>* segments) {
  ConcatenatedPins pins(fileNames);

  /* Create a temporary tab file to concatenate tab files in to */ 
  std::string inputFN = outputDir + "/converters-tmp.tcb";
  ofstream outFile;
  outFile.open(inputFN);
  if (!outFile.is_open()) {
    throw MyException("ERROR: Could not create the temporary file " + inputFN +
                      " for concatenating the input files.\n");
  }
  std::string line;
  while (pins.appendNextLine(line)) {
    outFile << line;
    line.clear();
  }
  outFile.close();
  if (!outFile) {
    throw MyException("ERROR: Could not write the temporary file " + inputFN +
                      ".\n");
  }
  if (segments) {
    *segments = pins.segments();
  }
  return inputFN;
}
