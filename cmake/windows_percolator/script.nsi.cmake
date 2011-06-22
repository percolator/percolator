!include "MUI2.nsh"
Name "percolator-@PACKAGE_VERSION@"
outfile "@CMAKE_BINARY_DIR@/percolator-@PACKAGE_VERSION@-win32.exe"
installDir $PROGRAMFILES\percolator

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_LICENSE "@CMAKE_SOURCE_DIR@/src/COPYING"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Section "Dummy Section" SecDummy
  SetOutPath "$INSTDIR"
  file @CMAKE_BINARY_DIR@/src/percolator.exe
# I had problems installing Tokyo Cabinet and XDR on Windows /Erik
#  file /tmp/percolatorBuild/src/converters/sqt2pin.exe
#  file @CMAKE_BINARY_DIR@/src/converters/mzidentml2pin.exe
#  file @CMAKE_BINARY_DIR@/src/elude/elude.exe
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Pin Xml Schema" SecPinXml
  SetOutPath "@PIN_SCHEMA_LOCATION@"
  file @CMAKE_SOURCE_DIR@/src/xml/percolator_in.xsd
SectionEnd

Section "Pout Xml Schema" SecPoutXml
  SetOutPath "@POUT_SCHEMA_LOCATION@"
  file @CMAKE_SOURCE_DIR@/src/xml/percolator_out.xsd
SectionEnd

#Section "Elude Models" SecElude
#  SetOutPath "@SCHEMA@..\elude_model"
#  file @CMAKE_SOURCE_DIR@/data/elude_model/
#SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\percolator.exe"
#  Delete "$INSTDIR\sqt2pin.exe"
#  Delete "$INSTDIR\mzidentml2pin.exe"
#  Delete "$INSTDIR\elude.exe"
  RMDir "$INSTDIR"
SectionEnd
