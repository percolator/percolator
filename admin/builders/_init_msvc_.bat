set BUILD_TARGET=%1

call :downloadfile https://github.com/microsoft/vswhere/releases/download/2.8.4/vswhere.exe vswhere.exe
for /f "usebackq tokens=*" %%i in (`.\vswhere.exe -legacy -latest  -property installationPath`) do (
  set MSVC_INSTALL_DIR=%%i
)

for /f "usebackq tokens=*" %%i in (`.\vswhere.exe -legacy -latest  -property installationVersion`) do (
  set MSVC_INSTALL_VERSION_FULL=%%i
)

set MSVC_VER=0
if %MSVC_INSTALL_VERSION_FULL% EQU 0 (
  echo Could not find a suitable Visual Studio version; supported versions: VS2013, VS2015, VS2017, VS2019
  EXIT /B 1  
) else (
  ::: strip off the main version number from the front of MSVC_INSTALL_VERSION_FULL
  set "MSVC_VER=%MSVC_INSTALL_VERSION_FULL:.=" & rem "%"
)

echo Visual Studio installation directory: %MSVC_INSTALL_DIR%
echo Visual Studio version: %MSVC_INSTALL_VERSION_FULL%
echo Visual Studio main version: %MSVC_VER%

set PROGRAM_FILES_DIR=C:\Program Files
set BUILD_PLATFORM=32bit
REG QUERY HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\VisualStudio\%MSVC_VER%.0\Setup > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo Build platform detected: 64-bit
  set BUILD_PLATFORM=64bit
  ::: double quotes around set command ensure that string is not evaluated
  set "PROGRAM_FILES_DIR=C:\Program Files (x86)"
)

::: The vcvarsall.bat has been moved in VS2017 and later
set AUXILIARY_PATH=
if %MSVC_VER% GEQ 15 (
  set AUXILIARY_PATH=Auxiliary\Build\
)

:: use the VS command prompt settings to set-up paths for compiler and builder
:: see https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx for possible vcvarsall.bat arguments
set VCVARS_BAT=%MSVC_INSTALL_DIR%\VC\%AUXILIARY_PATH%vcvarsall.bat
if not defined DevEnvDir (
  echo Setting compiler and builder paths for %BUILD_TARGET% target on %BUILD_PLATFORM% build platform
  if "%BUILD_TARGET%" == "32bit" (
    if "%BUILD_PLATFORM%" == "64bit" (
      call "%VCVARS_BAT%" amd64_x86
    ) else (
      call "%VCVARS_BAT%" x86
    )
  ) else (
    if "%BUILD_PLATFORM%" == "64bit" (
      call "%VCVARS_BAT%" amd64
    ) else (
      call "%VCVARS_BAT%" x86_amd64
    )
  )
)

EXIT /B %ERRORLEVEL%

:downloadfile
echo Downloading %1 to %2
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B
