@echo off
set MSVC_VER=0

:: use VS2015 if available
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.14.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo Using Visual Studio 2015
  set MSVC_VER=14
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

:: fall back to VS2013 is available
if %MSVC_VER% EQU 0 (
  REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.12.0 > nul 2> nul
  if %ERRORLEVEL% EQU 0 (
    echo Using Visual Studio 2013
    set MSVC_VER=12
  ) else (
    :: reset ERRORLEVEL to 0
    cd .
  )
)

if %MSVC_VER% EQU 0 (
  echo Could not find a suitable Visual Studio version; supported versions: VS2013, VS2015
  EXIT /B 1  
)

set PROGRAM_FILES_DIR=C:\Program Files
set BUILD_PLATFORM=32bit
REG QUERY HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\VisualStudio\%MSVC_VER%.0\Setup\VS > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo platform detected: 64-bit
  set BUILD_PLATFORM=64bit
  set "PROGRAM_FILES_DIR=C:\Program Files (x86)"
) else (
  :: reset ERRORLEVEL to 0
  cd .
)

set VCTARGET=%PROGRAM_FILES_DIR%\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build\win64
set RELEASE_DIR=%SRC_DIR%\release\win64
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
:: see https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx for possible vcvarsall.bat arguments
if not defined DevEnvDir (
  call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"
  if "%BUILD_PLATFORM%" == "64bit" (
    echo Setting variables for 64-bit
    call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" amd64
  ) else (
    echo Setting variables for 32-bit
    call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" x86_amd64
  )
)

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
  b2 address-model=64 threading=multi -j4 --with-system --with-filesystem --with-serialization -d0
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
CALL :getabspath PYTHON_DIR "%PYTHON_DIR%"
set PYTHON_URL=http://www.python.org/ftp/python/3.3.3/python-3.3.3.msi
if not exist "%PYTHON_DIR%" (
  echo Downloading and installing Python
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%PYTHON_URL%','%INSTALL_DIR%\python.msi')"
  cd /D "%INSTALL_DIR%"
  msiexec /i python.msi /quiet TARGETDIR="%PYTHON_DIR%" /Li python_install.log
)
setlocal
set PATH=%PATH%;%PYTHON_DIR%

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
set XERCES_DIR=%INSTALL_DIR%\xerces-c-3.1.1-x86_64-windows-vc-10.0
set XERCES_URL=http://archive.apache.org/dist/xerces/c/3/binaries/xerces-c-3.1.1-x86_64-windows-vc-10.0.zip
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
set SQLITE_DIR=%INSTALL_DIR%\sqlite3_x64
set SQLITE_SRC_URL=https://sqlite.org/2015/sqlite-amalgamation-3080803.zip
set SQLITE_DLL_URL=https://system.data.sqlite.org/blobs/1.0.96.0/sqlite-netFx45-binary-x64-2012-1.0.96.0.zip
if not exist "%SQLITE_DIR%" (
  echo Downloading and installing SQLite3
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_SRC_URL%','%INSTALL_DIR%\sqlite_src.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_DLL_URL%','%INSTALL_DIR%\sqlite_dll.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_src.zip" -o"%SQLITE_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_dll.zip" -o"%SQLITE_DIR%" > NUL
  
  ::: Generate lib from dll
  cd /D "%SQLITE_DIR%"
  ren sqlite-amalgamation-3080803 src
  ren SQLite.Interop.dll sqlite3.dll
  setlocal enableDelayedExpansion
  set DLL_BASE=%SQLITE_DIR%\sqlite3
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  lib /DEF:"!DEF_FILE!" /MACHINE:X64
  endlocal
)
set SQLITE_DIR=%SQLITE_DIR%;%SQLITE_DIR%\src

::: Needed for converters package and for system tests :::
set ZLIB_DIR=%INSTALL_DIR%\zlib_x64
set ZLIB_SRC_URL=http://win32builder.gnome.org/packages/3.6/zlib-dev_1.2.7-1_win64.zip
set ZLIB_DLL_URL=http://win32builder.gnome.org/packages/3.6/zlib_1.2.7-1_win64.zip
if not exist "%ZLIB_DIR%" (
  echo Downloading and installing ZLIB
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZLIB_SRC_URL%','%INSTALL_DIR%\zlib_src.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZLIB_DLL_URL%','%INSTALL_DIR%\zlib_dll.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\zlib_src.zip" -o"%ZLIB_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\zlib_dll.zip" -o"%ZLIB_DIR%" > NUL
  
  ::: Generate lib from dll
  setlocal enableDelayedExpansion
  set DLL_BASE=%ZLIB_DIR%\bin\zlib1
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  cd /D "%ZLIB_DIR%\bin"
  lib /DEF:"!DEF_FILE!" /MACHINE:X64
  endlocal
)
set ZLIB_DIR=%ZLIB_DIR%\bin;%ZLIB_DIR%\include
set PATH=%PATH%;%ZLIB_DIR%

::: needed for Elude :::
set DIRENT_H_PATH=%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\include\dirent.h
set DIRENT_H_URL=https://github.com/tronkko/dirent/archive/1.23.1.zip
if not exist "%DIRENT_H_PATH%" ( 
  echo Downloading and installing dirent.h 
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%DIRENT_H_URL%','%INSTALL_DIR%\dirent.zip')"
  %ZIP_EXE% x -aoa "%INSTALL_DIR%\dirent.zip" -o"%INSTALL_DIR%\dirent"
  copy "%INSTALL_DIR%\dirent\dirent-1.23.1\include\dirent.h" "%DIRENT_H_PATH%"
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
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DCMAKE_BUILD_TYPE=%BUILD_TYPE% -DXML_SUPPORT=OFF "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building percolator :::::::
if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")
cd /D "%BUILD_DIR%\percolator"
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DCMAKE_BUILD_TYPE=%BUILD_TYPE% -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building converters :::::::
if not exist "%BUILD_DIR%\converters" (md "%BUILD_DIR%\converters")
cd /D "%BUILD_DIR%\converters"
echo cmake converters.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DCMAKE_BUILD_TYPE=%BUILD_TYPE% -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%;%SQLITE_DIR%;%ZLIB_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator\src\converters"
echo build converters (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::: Building elude ::::::: 
if not exist "%BUILD_DIR%\elude" (md "%BUILD_DIR%\elude")
cd /D "%BUILD_DIR%\elude"
echo cmake elude.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DCMAKE_BUILD_TYPE=%BUILD_TYPE% -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" "%SRC_DIR%\percolator\src\elude_tool"
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

EXIT /B %errorlevel%

::: subroutines
:getabspath
SET "%1=%~f2"
EXIT /B
