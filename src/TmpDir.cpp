#include "TmpDir.h"

void TmpDir::createTempFile(std::string* tcf, char* tcd) {
  std::string str;

  try {
          boost::filesystem::path ph = boost::filesystem::unique_path();
          boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
          boost::filesystem::path file("converters-tmp.tcb");
          *tcf = std::string((dir / file).string());
          str =  dir.string();
          tcd = new char[str.size() + 1];
          std::copy(str.begin(), str.end(), tcd);
          tcd[str.size()] = '\0';
          if (boost::filesystem::is_directory(dir)) {
            boost::filesystem::remove_all(dir);
          }

          boost::filesystem::create_directory(dir);
        } catch (boost::filesystem::filesystem_error &e) {
          std::cerr << e.what() << std::endl;
        }
}