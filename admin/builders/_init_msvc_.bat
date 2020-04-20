set BUILD_TARGET=%1

set MSVC_VER=0

:: use VS2019 if available
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.16.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo Using Visual Studio 2019
  set MSVC_VER=16
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

:: use VS2017 if available (cannot check %ERRORLEVEL% inside an if statement!)
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.15.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  if %MSVC_VER% EQU 0 (
    echo Using Visual Studio 2017
    set MSVC_VER=15
  )
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

:: use VS2015 if available
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.14.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  if %MSVC_VER% EQU 0 (
    echo Using Visual Studio 2015
    set MSVC_VER=14
  )
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

:: use VS2013 if available
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.12.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  if %MSVC_VER% EQU 0 (
    echo Using Visual Studio 2013
    set MSVC_VER=12
  )
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

if %MSVC_VER% EQU 0 (
  echo Could not find a suitable Visual Studio version; supported versions: VS2013, VS2015, VS2017, VS2019
  EXIT /B 1
)

set PROGRAM_FILES_DIR=%ProgramFiles%
:::set PROGRAM_FILES_DIR=C:\Program Files
set BUILD_PLATFORM=32bit
::: REG QUERY HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\VisualStudio\%MSVC_VER%.0\Setup\VS > nul 2> nul
REG QUERY HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\VisualStudio\%MSVC_VER%.0\Setup > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo platform detected: 64-bit
  set BUILD_PLATFORM=64bit
  ::: double quotes around set command ensure that string is not evaluated
  set "PROGRAM_FILES_DIR=C:\Program Files (x86)"
) else (
  :: reset ERRORLEVEL to 0
  cd .
)

set VCTARGET=%PROGRAM_FILES_DIR%\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0

:: use the VS command prompt settings to set-up paths for compiler and builder
:: see https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx for possible vcvarsall.bat arguments
set VCVARS_BAT=%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat
if not defined DevEnvDir (
  call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"
  if "%BUILD_TARGET%" == "32bit" (
    if "%BUILD_PLATFORM%" == "64bit" (
      echo Setting variables for 64-bit
      call "%VCVARS_BAT%" amd64_x86
    ) else (
      echo Setting variables for 32-bit
      call "%VCVARS_BAT%" x86
    )
  ) else (
    if "%BUILD_PLATFORM%" == "64bit" (
      echo Setting variables for 64-bit
      call "%VCVARS_BAT%" amd64
    ) else (
      echo Setting variables for 32-bit
      call "%VCVARS_BAT%" x86_amd64
    )
  )
)
