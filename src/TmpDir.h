#ifndef TMPDIR_H_
#define TMPDIR_H_

#include <string>

/* A new directory in the system's temporary directory, removed with its
   content when the object is destroyed. */
class TmpDir {
  public:
    TmpDir();
    ~TmpDir();
    TmpDir(const TmpDir&) = delete;
    TmpDir& operator=(const TmpDir&) = delete;

    const std::string& getPath() const { return path_; }

  private:
    std::string path_;
};

#endif /* TMPDIR_H_ */
