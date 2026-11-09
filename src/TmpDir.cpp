#include "TmpDir.h"

#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>

#include "MyException.h"

TmpDir::TmpDir() {
  namespace fs = std::filesystem;
  try {
    std::random_device seed;
    std::mt19937_64 generator(seed());
    fs::path dir;
    do {
      std::ostringstream name;
      name << "percolator-" << std::hex << generator();
      dir = fs::temp_directory_path() / name.str();
    } while (!fs::create_directory(dir));
    path_ = dir.string();
  } catch (fs::filesystem_error& e) {
    throw MyException(std::string("ERROR: Could not create a temporary "
                                  "directory: ") + e.what() + "\n");
  }
}

TmpDir::~TmpDir() {
  std::error_code ec;  // never throw from a destructor
  std::filesystem::remove_all(path_, ec);
}
