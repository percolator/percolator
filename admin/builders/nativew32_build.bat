@echo off
set MSVC_VER=12
set VCTARGET=C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build-32bit
set RELEASE_DIR=%SRC_DIR%\release-32bit
set BUILD_TYPE=Release

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
SHIFT
GOTO parse
:endparse

:: use the VS command prompt settings to set-up paths for compiler and builder
call "C:\Program Files\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal
set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set ZIP_URL=http://downloads.sourceforge.net/sevenzip/7z920.exe
if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZIP_URL%','%INSTALL_DIR%\7zip.exe')"
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

set CMAKE_BASE=cmake-3.5.2-win32-x86
set CMAKE_URL=https://cmake.org/files/v3.5/%CMAKE_BASE%.zip
if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%CMAKE_URL%','%INSTALL_DIR%\cmake.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"

:: The windows binary release takes up 3GB, therefore we built only the libraries we need from source.
set BOOST_ROOT=%INSTALL_DIR%\boost_1_61_0
set BOOST_URL=http://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.7z/download
if not exist "%BOOST_ROOT%" (
  echo Downloading and installing Boost, this can take a few minutes...
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%BOOST_URL%','%INSTALL_DIR%\boost.7z')"
  %ZIP_EXE% x "%INSTALL_DIR%\boost.7z" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
  cd /D "%BOOST_ROOT%"
  call bootstrap
  bjam threading=multi -j4 --with-system --with-filesystem --with-serialization -d0
)
set BOOST_LIB=%BOOST_ROOT%\stage\lib

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS 3 Pre-release/3.0rc1/nsis-3.0rc1-setup.exe/download
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%NSIS_URL%','%INSTALL_DIR%\nsis.exe')"
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)
setlocal
set PATH=%PATH%;%INSTALL_DIR%\nsis

::: Needed for system tests :::
set PYTHON_DIR=%INSTALL_DIR%\python
set PYTHON_URL=http://www.python.org/ftp/python/3.3.3/python-3.3.3.msi
if not exist "%PYTHON_DIR%" (
  echo Downloading and installing Python
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%PYTHON_URL%','%INSTALL_DIR%\python.msi')"
  cd /D "%INSTALL_DIR%"
  msiexec /i python.msi /quiet TARGETDIR=python /Li python_install.log
)
setlocal
set PATH=%PATH%;%INSTALL_DIR%\python

::: Needed for system tests :::
set LIBXML_DIR=%INSTALL_DIR%\libxml2-2.7.8.win32
set LIBXML_URL=http://xmlsoft.org/sources/win32/libxml2-2.7.8.win32.zip
set ICONV_URL=http://sourceforge.net/projects/gettext/files/libiconv-win32/1.9.1/libiconv-1.9.1.bin.woe32.zip/download
set GETTEXT_URL=http://ftp.gnu.org/gnu/gettext/gettext-runtime-0.13.1.bin.woe32.zip
if not exist "%LIBXML_DIR%" (
  echo Downloading and installing LibXML
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%LIBXML_URL%','%INSTALL_DIR%\libxml.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ICONV_URL%','%INSTALL_DIR%\iconv.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%GETTEXT_URL%','%INSTALL_DIR%\gettext.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\libxml.zip" -o"%INSTALL_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\iconv.zip" -o"%LIBXML_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\gettext.zip" -o"%LIBXML_DIR%" > NUL
)
set PATH=%PATH%;%LIBXML_DIR%\bin

::: Needed for converters package and xml support in percolator package :::
set XERCES_DIR=%INSTALL_DIR%\xerces-c-3.1.1-x86-windows-vc-10.0
set XERCES_URL=http://archive.apache.org/dist/xerces/c/3/binaries/xerces-c-3.1.1-x86-windows-vc-10.0.zip
if not exist "%XERCES_DIR%" (
  echo Downloading and installing Xerces-C
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%XERCES_URL%','%INSTALL_DIR%\xerces.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\xerces.zip" -o"%INSTALL_DIR%" > NUL
)

::: Needed for converters package and xml support in percolator package :::
set XSD_DIR=%INSTALL_DIR%\xsd-3.3.0-i686-windows
set XSD_URL=http://www.codesynthesis.com/download/xsd/3.3/windows/i686/xsd-3.3.0-i686-windows.zip
if not exist "%XSD_DIR%" (
  echo Downloading and installing CodeSynthesis XSD
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%XSD_URL%','%INSTALL_DIR%\xsd.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\xsd.zip" -o"%INSTALL_DIR%" > NUL
)

::: Needed for converters package :::
set SQLITE_DIR=%INSTALL_DIR%\sqlite3
set SQLITE_SRC_URL=http://www.sqlite.org/2013/sqlite-amalgamation-3080200.zip
set SQLITE_DLL_URL=http://www.sqlite.org/2013/sqlite-dll-win32-x86-3080200.zip
if not exist "%SQLITE_DIR%" (
  echo Downloading and installing SQLite3
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_SRC_URL%','%INSTALL_DIR%\sqlite_src.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_DLL_URL%','%INSTALL_DIR%\sqlite_dll.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_src.zip" -o"%SQLITE_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_dll.zip" -o"%SQLITE_DIR%" > NUL
  cd /D "%SQLITE_DIR%"
  lib /DEF:"%SQLITE_DIR%\sqlite3.def" /MACHINE:X86
  ren sqlite-amalgamation-3080200 src
)
set SQLITE_DIR=%SQLITE_DIR%;%SQLITE_DIR%\src

::: Needed for converters package and for system tests :::
set ZLIB_DIR=%INSTALL_DIR%\zlib
set ZLIB_URL=http://sourceforge.net/projects/libpng/files/zlib/1.2.8/zlib128-dll.zip/download
if not exist "%ZLIB_DIR%" (
  echo Downloading and installing ZLIB
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZLIB_URL%','%INSTALL_DIR%\zlib.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\zlib.zip" -o"%ZLIB_DIR%" > NUL
  move "%ZLIB_DIR%\zlib1.dll" "%ZLIB_DIR%\lib\zlib1.dll" > NUL
)
set ZLIB_DIR=%ZLIB_DIR%;%ZLIB_DIR%\include
set PATH=%PATH%;%ZLIB_DIR%

::: needed for Elude :::
set DIRENT_H_PATH=C:\Program Files\Microsoft Visual Studio %MSVC_VER%.0\VC\include\dirent.h
set DIRENT_H_URL=http://www.softagalleria.net/download/dirent/dirent-1.20.1.zip
if not exist "%DIRENT_H_PATH%" ( 
  echo Downloading and installing dirent.h 
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%DIRENT_H_URL%','%INSTALL_DIR%\dirent.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\dirent.zip" -o"%INSTALL_DIR%\dirent" > NUL
  move "%INSTALL_DIR%\dirent\include\dirent.h" "%DIRENT_H_PATH%" > NUL
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: END INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::
:::::::::::: START BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::

if not exist "%BUILD_DIR%" (md "%BUILD_DIR%")

::::::: Building percolator without xml support :::::::
if not exist "%BUILD_DIR%\percolator-noxml" (md "%BUILD_DIR%\percolator-noxml")
cd /D "%BUILD_DIR%\percolator-noxml"
echo cmake percolator-noxml.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DXML_SUPPORT=OFF "%SRC_DIR%\percolator"
echo build percolator-noxml (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building percolator :::::::
if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")
cd /D "%BUILD_DIR%\percolator"
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building converters :::::::
if not exist "%BUILD_DIR%\converters" (md "%BUILD_DIR%\converters")
cd /D "%BUILD_DIR%\converters"
echo cmake converters.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%;%SQLITE_DIR%;%ZLIB_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator\src\converters"
echo build converters (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building elude ::::::: 
if not exist "%BUILD_DIR%\elude" (md "%BUILD_DIR%\elude")
cd /D "%BUILD_DIR%\elude"
echo cmake elude.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" "%SRC_DIR%\percolator\src\elude_tool"
echo build elude (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
copy "%BUILD_DIR%\percolator-noxml\per*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\percolator\per*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\converters\per*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\elude\elude*.exe" "%RELEASE_DIR%"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%\percolator\admin\builders"
