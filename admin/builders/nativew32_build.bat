@echo off
set MSVC_VER=12
set BOOST_ROOT=C:\local\boost_1_55_0
set BOOST_LIB=%BOOST_ROOT%\stage\lib
set XERCES_DIR=C:\local\xerces-c-3.1.1-x86-windows-vc-10.0
set XSD_DIR=C:\local\xsd-3.3.0-i686-windows
set SQLITE_DIR=C:\local\sqlite3;C:\local\sqlite3-src
set ZLIB_DIR=C:\local\zlib128;C:\local\zlib128\include
set VCTARGET=C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=Z:\percolator
set BUILD_DIR=Z:\build

:: use the VS command prompt settings to set-up paths for compiler and builder
call "C:\Program Files\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"

md "%BUILD_DIR%"

md "%BUILD_DIR%\percolator"
cd "%BUILD_DIR%\percolator"

::::::: Building percolator :::::::
echo "cmake percolator....."
cmake -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%"
echo "build percolator (this will take a few minutes)....."
::msbuild PERCOLATOR.sln /p:VCTargetsPath="%VCTARGET%" /p:Configuration=Release
msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=Release /m

md "%BUILD_DIR%\converters"
cd "%BUILD_DIR%\converters"

::::::: Building percolator :::::::
echo "cmake converters....."
cmake -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%;%SQLITE_DIR%;%ZLIB_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\src\converters"
echo "build converters (this will take a few minutes)....."
::msbuild PERCOLATOR.sln /p:VCTargetsPath="%VCTARGET%" /p:Configuration=Release
msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=Release /m


cd %SRC_DIR%\admin\builders
