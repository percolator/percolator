!include "MUI2.nsh"
Name "elude-@PACKAGE_VERSION@"
outfile "@CMAKE_BINARY_DIR@/elude-@PACKAGE_VERSION@-win32.exe"
installDir $PROGRAMFILES\elude

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_LICENSE "@CMAKE_SOURCE_DIR@/COPYING"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Section "Dummy Section" SecDummy
  SetOutPath "$INSTDIR"
  file @CMAKE_BINARY_DIR@/elude.exe
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Elude Models" SecElude
  SetOutPath "@ELUDE_MODELS_PATH@elude_model"
  file @CMAKE_SOURCE_DIR@/data/elude/models/
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\elude.exe"
  RMDir "$INSTDIR"
SectionEnd
