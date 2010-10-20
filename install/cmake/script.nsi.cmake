!include "MUI2.nsh"
Name "percolator-@PACKAGE_VERSION@"
outfile "@CMAKE_BINARY_DIR@/percolator-@PACKAGE_VERSION@-win32.exe"
installDir $PROGRAMFILES\percolator

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_LICENSE "@CMAKE_SOURCE_DIR@/install/cmake/license-win32-statically-built"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Section "Dummy Section" SecDummy
  SetOutPath "$INSTDIR"
  file tmp/percolatorInstall/percolator.exe
# I had problems installing Tokyo Cabinet and XDR on Windows /Erik
#  file src/converters/sqt2pin.exe
#  file src/converters/mzidentml2pin.exe
#  file src/elude/elude.exe
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\percolator.exe"
#  Delete "$INSTDIR\sqt2pin.exe"
#  Delete "$INSTDIR\mzidentml2pin.exe"
#  Delete "$INSTDIR\elude.exe"
  RMDir "$INSTDIR"
SectionEnd
