@echo off
:: Builds Percolator for 32-bit Windows with Visual Studio and makes the NSIS
:: installer. Three cmake commands; nothing is downloaded or installed.
::
:: Needs: CMake 3.15 or later on the PATH, Visual Studio 2019 or later with
:: the C++ workload (CMake finds the newest one), and NSIS for the installer.
:: The GitHub runner images hold all three.
::
:: usage: nativew32_build.bat [-s <dir that holds percolator\>] [-b <build dir>] [-r <release dir>]
::   The build tree is <build dir>\percolator; the installer is copied to
::   <release dir>.

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

where cmake > nul 2> nul
if errorlevel 1 (
  echo ERROR: cmake is not on the PATH
  exit /B 1
)

if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")
if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")

echo cmake percolator.....
cmake -S "%SRC_DIR%\percolator" -B "%BUILD_DIR%\percolator" -A Win32
if errorlevel 1 exit /B 1

echo build percolator (this will take a few minutes).....
cmake --build "%BUILD_DIR%\percolator" --config %BUILD_TYPE% --parallel
if errorlevel 1 exit /B 1

echo package percolator.....
cmake --build "%BUILD_DIR%\percolator" --config %BUILD_TYPE% --target PACKAGE
if errorlevel 1 exit /B 1

echo Copying installers to %RELEASE_DIR%
copy /Y "%BUILD_DIR%\percolator\per*.exe" "%RELEASE_DIR%"
if errorlevel 1 exit /B 1

echo Finished buildscript execution in build directory %BUILD_DIR%
exit /B 0
