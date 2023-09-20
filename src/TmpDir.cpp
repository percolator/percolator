#include "TmpDir.h"

void TmpDir::createTempFile(std::string& tcf, std::string& tcd) {
  //TODO it would be nice to somehow avoid these declararions and therefore avoid the linking to
  //boost filesystem when we don't use them      
  try {
          boost::filesystem::path ph = boost::filesystem::unique_path();
          boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
          boost::filesystem::path file("converters-tmp.tcb");
          tcf = std::string((dir / file).string());
          tcd = dir.string();
          if (boost::filesystem::is_directory(dir)) {
            boost::filesystem::remove_all(dir);
          }

          boost::filesystem::create_directory(dir);
        } catch (boost::filesystem::filesystem_error &e) {
          std::cerr << e.what() << std::endl;
        }
}

