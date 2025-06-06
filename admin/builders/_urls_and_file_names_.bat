::: Centralized place for urls and files for all windows builders ...
::: please do not change compression type in urls, since decompression is
::: hardcoded in the respective buiding scripts

::: 7-zip
set ZIP_URL=https://www.7-zip.org/a/7z2301.exe

::: CMake
set CMAKE_VERSION=3.31.6
set CMAKE_BASE=cmake-%CMAKE_VERSION%-windows-x86_64
set CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v%CMAKE_VERSION%/%CMAKE_BASE%.zip

::: Boost
set BOOST_BASE=boost_1_71_0
set BOOST_URL=https://archives.boost.io/release/1.71.0/source/boost_1_71_0.7z

::: NSIS
set NSIS_URL=https://prdownloads.sourceforge.net/nsis/nsis-3.11-setup.exe?download

::: Python
set PYTHON_URL=https://www.python.org/ftp/python/3.4.4/python-3.4.4.msi

::: Libxml
set LIBXML_BASE=libxml2-2.7.8.win32
set LIBXML_URL=http://xmlsoft.org/sources/win32/%LIBXML_BASE%.zip
set ICONV_URL=https://downloads.sourceforge.net/project/gettext/libiconv-win32/1.9.1/libiconv-1.9.1.bin.woe32.zip
set GETTEXT_URL=https://ftp.gnu.org/gnu/gettext/gettext-runtime-0.13.1.bin.woe32.zip

::: XercesC
set XERCES_BASE=xerces-c-3.1.1-x86-windows-vc-10.0
set XERCES_URL=https://archive.apache.org/dist/xerces/c/3/binaries/%XERCES_BASE%.zip
set XERCES_64_BASE=xerces-c-3.1.1-x86_64-windows-vc-10.0
set XERCES_64_URL=https://archive.apache.org/dist/xerces/c/3/binaries/%XERCES_64_BASE%.zip

::: XSD
::: set XSD_BASE=xsd-3.3.0-i686-windows
::: set XSD_URL=https://www.codesynthesis.com/download/xsd/3.3/windows/i686/%XSD_BASE%.zip
set XSD_BASE=xsd-4.2.0-x86_64-windows10
set XSD_URL=https://www.codesynthesis.com/download/xsd/4.2/windows/windows10/x86_64/%XSD_BASE%.zip
set LIBXSD_BASE=libxsd-4.2.0-windows
set LIBXSD_URL=https://www.codesynthesis.com/download/xsd/4.2/windows/windows10/x86_64/%LIBXSD_BASE%.zip

::: SQLite3
set SQLITE_SRC_BASE=sqlite-amalgamation-3080200
set SQLITE_SRC_URL=https://www.sqlite.org/2013/%SQLITE_SRC_BASE%.zip
set SQLITE_DLL_URL=https://www.sqlite.org/2013/sqlite-dll-win32-x86-3080200.zip

set SQLITE_64_SRC_BASE=sqlite-amalgamation-3080803
set SQLITE_64_SRC_URL=https://sqlite.org/2015/%SQLITE_64_SRC_BASE%.zip
set SQLITE_64_DLL_URL=https://system.data.sqlite.org/downloads/1.0.119.0/sqlite-netFx45-binary-x64-2012-1.0.119.0.zip
::: Zlib
set ZLIB_URL=https://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib128-dll.zip

set ZLIB_64_URL=https://nsis.sourceforge.io/mediawiki/images/b/bb/Zlib-1.2.8-win64-AMD64.zip

::: dirent.h
set DIRENT_H_VERSION=1.23.1
set DIRENT_H_URL=https://github.com/tronkko/dirent/archive/%DIRENT_H_VERSION%.zip
