net use z: \\vboxsrv\vagrant

@echo off

setlocal

set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build\win32
set RELEASE_DIR=%SRC_DIR%\release\win32
set BUILD_TYPE=Release

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
SHIFT
GOTO parse
:endparse

cd /D "%SRC_DIR%"

call percolator\admin\builders\_init_msvc_.bat 32bit
if %errorlevel% NEQ 0 (
  EXIT /B %errorlevel%
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

call percolator\admin\builders\_urls_and_file_names_.bat

set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  call :downloadfile %ZIP_URL% %INSTALL_DIR%\7zip.exe
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"

:: The windows binary release takes up 3GB, therefore we build only the libraries we need from source.
set BOOST_ROOT=%INSTALL_DIR%\%BOOST_BASE%
if not exist "%BOOST_ROOT%" (
  echo Downloading and installing Boost, this can take a few minutes...
  call :downloadfile %BOOST_URL% %INSTALL_DIR%\boost.7z
  %ZIP_EXE% x "%INSTALL_DIR%\boost.7z" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
  cd /D "%BOOST_ROOT%"
  call bootstrap
  b2 threading=multi -j4 --with-system --with-filesystem --with-serialization -d0
)
set BOOST_LIB=%BOOST_ROOT%\stage\lib

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  call :downloadfile "%NSIS_URL%" %INSTALL_DIR%\nsis.exe
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)
setlocal
set PATH=%PATH%;%INSTALL_DIR%\nsis

::: Needed for system tests :::
set PYTHON_DIR=%INSTALL_DIR%\python
CALL :getabspath PYTHON_DIR "%PYTHON_DIR%"
if not exist "%PYTHON_DIR%" (
  echo Downloading and installing Python
  call :downloadfile %PYTHON_URL% %INSTALL_DIR%\python.msi
  cd /D "%INSTALL_DIR%"
  msiexec /i python.msi /quiet TARGETDIR="%PYTHON_DIR%" /Li python_install.log
)
setlocal
set PATH=%PATH%;%INSTALL_DIR%\python

::: Needed for system tests :::
set LIBXML_DIR=%INSTALL_DIR%\%LIB_XML_BASE%
if not exist "%LIBXML_DIR%" (
  echo Downloading and installing LibXML
  call :downloadfile %LIBXML_URL% %INSTALL_DIR%\libxml.zip
  call :downloadfile %ICONV_URL% %INSTALL_DIR%\iconv.zip
  call :downloadfile %GETTEXT_URL% %INSTALL_DIR%\gettext.zip
  %ZIP_EXE% x "%INSTALL_DIR%\libxml.zip" -o"%INSTALL_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\iconv.zip" -o"%LIBXML_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\gettext.zip" -o"%LIBXML_DIR%" > NUL
)
set PATH=%PATH%;%LIBXML_DIR%\bin

::: Needed for converters package and xml support in percolator package :::
set XERCES_DIR=%INSTALL_DIR%\%XERCES_BASE%
if not exist "%XERCES_DIR%" (
  echo Downloading and installing Xerces-C
  call :downloadfile %XERCES_URL% %INSTALL_DIR%\xerces.zip
  %ZIP_EXE% x "%INSTALL_DIR%\xerces.zip" -o"%INSTALL_DIR%" > NUL
)

::: Needed for converters package and xml support in percolator package :::
set XSD_DIR=%INSTALL_DIR%\%XSD_BASE%
if not exist "%XSD_DIR%" (
  echo Downloading and installing CodeSynthesis XSD
  call :downloadfile %XSD_URL% %INSTALL_DIR%\xsd.zip
  %ZIP_EXE% x "%INSTALL_DIR%\xsd.zip" -o"%INSTALL_DIR%" > NUL
)

::: Needed for converters package :::
set SQLITE_DIR=%INSTALL_DIR%\sqlite3
if not exist "%SQLITE_DIR%" (
  echo Downloading and installing SQLite3
  call :downloadfile %SQLITE_SRC_URL% %INSTALL_DIR%\sqlite_src.zip
  call :downloadfile %SQLITE_DLL_URL% %INSTALL_DIR%\sqlite_dll.zip
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_src.zip" -o"%SQLITE_DIR%" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_dll.zip" -o"%SQLITE_DIR%" > NUL
  cd /D "%SQLITE_DIR%"
  lib /DEF:"%SQLITE_DIR%\sqlite3.def" /MACHINE:X86
  ren %SQLITE_SRC_BASE% src
)
set SQLITE_DIR=%SQLITE_DIR%;%SQLITE_DIR%\src

::: Needed for converters package and for system tests :::
set ZLIB_DIR=%INSTALL_DIR%\zlib
if not exist "%ZLIB_DIR%" (
  echo Downloading and installing ZLIB
  call :downloadfile %ZLIB_URL% %INSTALL_DIR%\zlib.zip
  %ZIP_EXE% x "%INSTALL_DIR%\zlib.zip" -o"%ZLIB_DIR%" > NUL
  move "%ZLIB_DIR%\zlib1.dll" "%ZLIB_DIR%\lib\zlib1.dll" > NUL
)
set ZLIB_DIR=%ZLIB_DIR%;%ZLIB_DIR%\include
set PATH=%PATH%;%ZLIB_DIR%

::: needed for  :::
set DIRENT_H_PATH=%VCToolsInstallDir%\include\dirent.h
if not exist "%DIRENT_H_PATH%" (
  echo Downloading and installing dirent.h
  call :downloadfile %DIRENT_H_URL% %INSTALL_DIR%\dirent.zip
  %ZIP_EXE% x -aoa "%INSTALL_DIR%\dirent.zip" -o"%INSTALL_DIR%" > NUL
  copy "%INSTALL_DIR%\dirent-%DIRENT_H_VERSION%\include\dirent.h" "%DIRENT_H_PATH%" > NUL
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
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m

::::::: Building percolator :::::::
if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")
cd /D "%BUILD_DIR%\percolator"
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::::::: Building converters :::::::
if not exist "%BUILD_DIR%\converters" (md "%BUILD_DIR%\converters")
cd /D "%BUILD_DIR%\converters"
echo cmake converters.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%;%SQLITE_DIR%;%ZLIB_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator\src\converters"
echo build converters (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
set /A exit_code=0
call :copytorelease "%BUILD_DIR%\percolator-noxml\per*.exe"
call :copytorelease "%BUILD_DIR%\percolator\per*.exe"
call :copytorelease "%BUILD_DIR%\converters\per*.exe"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%\percolator\admin\builders"

EXIT /B %exit_code%

::: subroutines
:getabspath
SET "%1=%~f2"
EXIT /B

:downloadfile
echo Downloading "%1" to "%2"
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B

:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
