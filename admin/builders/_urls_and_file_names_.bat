::: Centralized place for urls and files for all windows builders ...
::: please do not change compression type in urls, since decompression is
::: hardcoded in the respective buiding scripts

::: 7-zip
set ZIP_URL=https://www.7-zip.org/a/7z2301.exe

::: CMake
set CMAKE_VERSION=3.31.6
set CMAKE_BASE=cmake-%CMAKE_VERSION%-windows-x86_64
set CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v%CMAKE_VERSION%/%CMAKE_BASE%.zip

::: NSIS
set NSIS_URL=https://prdownloads.sourceforge.net/nsis/nsis-3.11-setup.exe?download

::: Python
set PYTHON_URL=https://www.python.org/ftp/python/3.4.4/python-3.4.4.msi

::: Zlib
set ZLIB_URL=https://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib128-dll.zip

set ZLIB_64_URL=https://nsis.sourceforge.io/mediawiki/images/b/bb/Zlib-1.2.8-win64-AMD64.zip

::: dirent.h
set DIRENT_H_VERSION=1.23.1
set DIRENT_H_URL=https://github.com/tronkko/dirent/archive/%DIRENT_H_VERSION%.zip
