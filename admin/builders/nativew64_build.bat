@echo off
set MSVC_VER=12
set VCTARGET=C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build
set RELEASE_DIR=%SRC_DIR%
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

set CMAKE_URL=http://www.cmake.org/files/v2.8/cmake-2.8.12.1-win32-x86.exe
if not exist "%INSTALL_DIR%\cmake" (
  echo Downloading and installing CMake
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%CMAKE_URL%','%INSTALL_DIR%\cmake.exe')"
  "%INSTALL_DIR%\cmake.exe" /S /D=%INSTALL_DIR%\cmake
)
set CMAKE_EXE="%INSTALL_DIR%\cmake\bin\cmake.exe"

set ZIP_URL=http://downloads.sourceforge.net/sevenzip/7z920.exe
if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZIP_URL%','%INSTALL_DIR%\7zip.exe')"
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

:: The windows binary release of boost 1.55 does not include the serialization library, therefore we built from source.
:: If this bug is fixed in the next version, we can just grab the binaries (much faster, except for a bigger download)
set BOOST_ROOT=%INSTALL_DIR%\boost_1_55_0
set BOOST_URL=http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.7z/download
if not exist "%BOOST_ROOT%" (
  echo Downloading and installing Boost, this can take a few minutes...
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%BOOST_URL%','%INSTALL_DIR%\boost.7z')"
  %ZIP_EXE% x "%INSTALL_DIR%\boost.7z" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
  :::::: Fix a bug in the 1.55 release ::::::::
  echo #include ^<algorithm^> > tmp.txt
  type "%BOOST_ROOT%\libs\serialization\src\basic_text_iprimitive.cpp" >> tmp.txt
  type tmp.txt > "%BOOST_ROOT%\libs\serialization\src\basic_text_iprimitive.cpp"
  echo #include ^<algorithm^> > tmp.txt
  type "%BOOST_ROOT%\libs\serialization\src\basic_text_oprimitive.cpp" >> tmp.txt
  type tmp.txt > "%BOOST_ROOT%\libs\serialization\src\basic_text_oprimitive.cpp"
  echo #include ^<algorithm^> > tmp.txt
  type "%BOOST_ROOT%\libs\serialization\src\basic_text_wiprimitive.cpp" >> tmp.txt
  type tmp.txt > "%BOOST_ROOT%\libs\serialization\src\basic_text_wiprimitive.cpp"
  echo #include ^<algorithm^> > tmp.txt
  type "%BOOST_ROOT%\libs\serialization\src\basic_text_woprimitive.cpp" >> tmp.txt
  type tmp.txt > "%BOOST_ROOT%\libs\serialization\src\basic_text_woprimitive.cpp"
  del tmp.txt
  :::::: end bug fix ::::::::
  cd /D "%BOOST_ROOT%"
  call bootstrap
  bjam address-model=64 threading=multi --with-system --with-filesystem --with-serialization -d0
)
set BOOST_LIB=%BOOST_ROOT%\stage\lib

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=http://downloads.sourceforge.net/project/nsis/NSIS 3 Pre-release/3.0a2/nsis-3.0a2-setup.exe
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%NSIS_URL%','%INSTALL_DIR%\nsis.exe')"
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)

::: Needed for system tests :::
set PYTHON_DIR=%INSTALL_DIR%\python
set PYTHON_URL=http://www.python.org/ftp/python/3.3.3/python-3.3.3.msi
if not exist "%PYTHON_DIR%" (
  echo Downloading and installing Python
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%PYTHON_URL%','%INSTALL_DIR%\python.msi')"
  msiexec /i "%INSTALL_DIR%\python.msi" /quiet TARGETDIR="%PYTHON_DIR%" /li "%INSTALL_DIR%\python_install.log"
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
set XERCES_URL=http://apache.mirrors.spacedump.net/xerces/c/3/binaries/xerces-c-3.1.1-x86_64-windows-vc-10.0.zip
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
set SQLITE_SRC_URL=http://www.sqlite.org/snapshot/sqlite-amalgamation32k-201409200035.zip
set SQLITE_DLL_URL=http://www.sqlite.org/snapshot/sqlite-dll-win64-x64-201409200035.zip
if not exist "%SQLITE_DIR%" (
  echo Downloading and installing SQLite3
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_SRC_URL%','%INSTALL_DIR%\sqlite_src.zip')"
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SQLITE_DLL_URL%','%INSTALL_DIR%\sqlite_dll.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_src.zip" -o"%SQLITE_DIR%\src" > NUL
  %ZIP_EXE% x "%INSTALL_DIR%\sqlite_dll.zip" -o"%SQLITE_DIR%" > NUL
  
  ::: Generate lib from dll
  setlocal enableDelayedExpansion
  set DLL_BASE=%SQLITE_DIR%\sqlite3
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  cd /D "%SQLITE_DIR%"
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

::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: END INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::
:::::::::::: START BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::

if not exist "%BUILD_DIR%" (md "%BUILD_DIR%")

if not exist "%BUILD_DIR%\percolator-noxml" (md "%BUILD_DIR%\percolator-noxml")
cd /D "%BUILD_DIR%\percolator-noxml"

::::::: Building percolator without xml support :::::::
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=OFF "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

if not exist "%BUILD_DIR%\percolator" (md "%BUILD_DIR%\percolator")
cd /D "%BUILD_DIR%\percolator"

::::::: Building percolator :::::::
echo cmake percolator.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator"
echo build percolator (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

if not exist "%BUILD_DIR%\converters" (md "%BUILD_DIR%\converters")
cd /D "%BUILD_DIR%\converters"

::::::: Building converters :::::::
echo cmake converters.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="%XERCES_DIR%;%XSD_DIR%;%SQLITE_DIR%;%ZLIB_DIR%" -DXML_SUPPORT=ON "%SRC_DIR%\percolator\src\converters"
echo build converters (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::if not exist "%BUILD_DIR%\elude" (md "%BUILD_DIR%\elude")
::cd /D "%BUILD_DIR%\elude"

::::::: Building elude (Not working at the moment, see https://github.com/percolator/percolator/issues/106)::::::: 
::echo cmake elude.....
::%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%BOOST_ROOT%" -DBOOST_LIBRARYDIR="%BOOST_LIB%" "%SRC_DIR%\percolator\src\elude_tool"
::echo build elude (this will take a few minutes).....
::msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

copy "%BUILD_DIR%\percolator-noxml\per*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\percolator\per*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\converters\per*.exe" "%RELEASE_DIR%"
::copy "%BUILD_DIR%\elude\elude*.exe" "%RELEASE_DIR%"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%\percolator\admin\builders"
