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

/*
 * Several tab-delimited input files as one input: the header of the file
 * with the most columns, then the PSM lines of the files in turn, with zeros
 * in the columns that a file lacks and without the DefaultDirection line of
 * any file but the first. The lines are those of the file that
 * ValidateTabFile::concatenateMultiplePINs writes, made one at a time
 * instead of copied into a temporary file of the size of the input.
 */
class ConcatenatedPins {
 public:
  // Reads the header and the first rows of every file: see the definition
  // for the messages and the error.
  explicit ConcatenatedPins(const std::vector<std::string>& fileNames);
  // Appends the next line, with its line break, to out. Returns false after
  // the last line.
  bool appendNextLine(std::string& out);
  // starts again at the first line
  void rewind();
  // One entry per input file with PSM lines whose lines have been reached,
  // in the order of the input: complete after one pass. A reader of the
  // lines can name the file and line of every line it has got.
  const std::vector<InputSegment>& segments() const { return segments_; }

 private:
  struct File {
    std::string name;
    std::vector<int> missingCols;  // columns of the header that it lacks
    unsigned int firstRowLine = 2u;  // the line of its first PSM
  };
  bool openNextFile();

  std::string headerLine_;
  std::vector<File> files_;  // the files with PSM lines
  std::vector<InputSegment> segments_;
  std::ifstream current_;
  std::string row_;
  std::size_t nextFile_ = 0u;
  unsigned int outLines_ = 0u;  // lines given so far
  bool headerDone_ = false, fileOpen_ = false;
};

/*
 * The lines of ConcatenatedPins as an input stream. seekg(0) starts again
 * at the first line; no other position can be sought.
 */
class ConcatenatedPinStream : public std::istream {
 public:
  explicit ConcatenatedPinStream(const std::vector<std::string>& fileNames)
      : std::istream(NULL), pins_(fileNames), buffer_(pins_) {
    rdbuf(&buffer_);
  }
  const std::vector<InputSegment>& segments() const {
    return pins_.segments();
  }

 private:
  class Buffer : public std::streambuf {
   public:
    explicit Buffer(ConcatenatedPins& pins) : pins_(pins) {}

   protected:
    int_type underflow() override;
    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                     std::ios_base::openmode which) override;
    pos_type seekpos(pos_type pos, std::ios_base::openmode which) override;

   private:
    static constexpr std::size_t kBlockSize = 1u << 16;
    ConcatenatedPins& pins_;
    std::string text_;
  };
  ConcatenatedPins pins_;
  Buffer buffer_;
};

class ValidateTabFile {
  public:
    // Writes the lines of ConcatenatedPins to a file in outputDir and
    // returns its name. If segments is given, it receives one entry per
    // input file with PSM lines, in the order of the concatenated file.
    std::string concatenateMultiplePINs(const std::vector<std::string>& fileNames,
                                        const std::string& outputDir,
                                        std::vector<InputSegment>* segments = nullptr);
};

#endif /*VALIDATETABFILE_H_*/
