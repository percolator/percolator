@echo off

setlocal

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

cd /D "%SRC_DIR%"

call percolator\admin\builders\_init_msvc_.bat 64bit
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

set ZIP_INSTALLER=%INSTALL_DIR%\7zip-installer.exe
set ZIP_EXE=%INSTALL_DIR%\7zip\7z.exe

if not exist "%ZIP_EXE%" (
  echo Downloading and installing 7-Zip...
  call :downloadfile %ZIP_URL% %ZIP_INSTALLER%

  :: Run installer silently to extract to the desired location
  start /wait "" "%ZIP_INSTALLER%" /S /D=%INSTALL_DIR%\7zip

  :: Verify the install
  if not exist "%ZIP_EXE%" (
    echo ERROR: 7z.exe not found after installation.
    exit /B 1
  )
)

if not exist "%ZIP_EXE%" (
  echo ERROR: 7z.exe not found after unpacking.
  exit /B 1
)

if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc -bso0 || (
    echo Extraction failed for cmake.zip
    EXIT /B 1
  )
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"





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
set PATH=%PATH%;%PYTHON_DIR%

::: Needed for converters package and for system tests :::
set ZLIB_DIR=%INSTALL_DIR%\zlib_x64
if not exist "%ZLIB_DIR%" (
  echo Downloading and installing ZLIB
  call :downloadfile %ZLIB_64_URL% %INSTALL_DIR%\zlib.zip
  %ZIP_EXE% x "%INSTALL_DIR%\zlib.zip" -o"%ZLIB_DIR%" > NUL
  mkdir "%ZLIB_DIR%\bin" > NUL
  move "%ZLIB_DIR%\zlib1.dll" "%ZLIB_DIR%\bin\zlib.dll" > NUL
)
set ZLIB_DIR=%ZLIB_DIR%\lib;%ZLIB_DIR%\include;%ZLIB_DIR%\bin
set PATH=%PATH%;%ZLIB_DIR%

::: needed for Elude :::
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

::::::: Building percolator :::::::
if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")
cd /D "%BUILD_DIR%\percolator"
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -A x64 -DCMAKE_BUILD_TYPE=%BUILD_TYPE% -DGOOGLE_TEST=1 ^
  "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
set /A exit_code=0
call :copytorelease "%BUILD_DIR%\percolator\per*.exe"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%\percolator\admin\builders"

EXIT /B %exit_code%

::: subroutines
:getabspath
SET "%1=%~f2"
EXIT /B

:downloadfile
echo Downloading "%1" to "%2"
curl -fsSL -o "%2" "%1"
if errorlevel 1 (
  echo ERROR: Failed to download %1
  exit /b 1
)
exit /b

:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
