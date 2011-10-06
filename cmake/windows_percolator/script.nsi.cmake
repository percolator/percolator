; Percolator installer


;-----------------------------------------------------------------------------
; Some paths.
;-----------------------------------------------------------------------------

!ifndef MING_PATH
    !define MING_PATH "/usr/i686-pc-mingw32/sys-root/mingw"
!endif
!define MING_BIN "${MING_PATH}/bin"
!define MING_LIB "${MING_PATH}/lib"
!define BUILD_PATH "@CMAKE_BINARY_DIR@"
!define SOURCE_PATH "@CMAKE_SOURCE_DIR@"
!define NSI_PATH "${SOURCE_PATH}/admin/win/nsi"


; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "Percolator"
!define PRODUCT_VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"
!define PRODUCT_PUBLISHER "<Lukas Kall>"
!define PRODUCT_WEB_SITE "http://per-colator.com/"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
!define PRODUCT_UNINST_ROOT_KEY "HKLM"
!define EXECUTABLE "percolator.exe"
!define PROGICON "${NSI_PATH}\installer.ico"
!define SETUP_BITMAP "${NSI_PATH}\welcome.bmp"


; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON ${PROGICON}
!define MUI_UNICON ${PROGICON}  ; define uninstall icon to appear in Add/Remove Programs

; Welcome page
!insertmacro MUI_PAGE_WELCOME
# ; License page
# !insertmacro MUI_PAGE_LICENSE "license.txt"   ; text file with license terms
; Directory page
!insertmacro MUI_PAGE_DIRECTORY
; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\readme.txt"  ; readme.txt file for user
!insertmacro MUI_PAGE_FINISH

; Uninstaller pages
!insertmacro MUI_UNPAGE_INSTFILES

; Language files
!insertmacro MUI_LANGUAGE "English"
; MUI end ------



;-----------------------------------------------------------------------------
; Installer version
;-----------------------------------------------------------------------------

!define VER_MAJOR "@CPACK_PACKAGE_VERSION_MAJOR@"
!define VER_MINOR "@CPACK_PACKAGE_VERSION_MINOR@"
!define VER_BUILD "@CPACK_PACKAGE_VERSION_PATCH@"
!define VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"
!define REVISION ""

Name "percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@"
Caption "Percolator Installer"
outfile "@CMAKE_BINARY_DIR@/percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@-win32.exe"
installDir $PROGRAMFILES\Percolator
InstallDirRegKey HKCU "Software\Percolator" ""
ShowInstDetails show
ShowUnInstDetails show



Section "MainSection" SEC01
   SectionIn 1 2 3 RO
   SetDetailsPrint listonly

   SetDetailsPrint textonly
   DetailPrint "Installing Percolator essentials."
   SetDetailsPrint listonly
   SetOutPath "$INSTDIR"
   File "${BUILD_PATH}\src\percolator.exe"
   
   ;MINGW
   File "${MING_BIN}\boost_filesystem-gcc45-mt-1_46_1.dll"
   File "${MING_BIN}\boost_system-gcc45-mt-1_46_1.dll"
   File "${MING_BIN}\libgcc_s_sjlj-1.dll"
   File "${MING_BIN}\libstdc++-6.dll"
   File "${MING_BIN}\libxerces-c-3-0.dll"

   ;License & release notes.
   File "@CPACK_RESOURCE_FILE_LICENSE@"
   File /oname=NOTES.txt ${NSI_PATH}\RELEASE_NOTES.txt

  CreateDirectory "$SMPROGRAMS\${PRODUCT_NAME}"
  SetShellVarContext all
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\${PRODUCT_NAME}.lnk" "$INSTDIR\${EXECUTABLE}" "" "$INSTDIR\${PROGICON}" 0
  CreateShortCut "$DESKTOP\${PRODUCT_NAME}.lnk" "$INSTDIR\${EXECUTABLE}" "" "$INSTDIR\${PROGICON}" 0
SectionEnd

Section -AdditionalIcons
  WriteIniStr "$INSTDIR\${PRODUCT_NAME}.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Website.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  SetShellVarContext all		; scope is "All Users"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall.lnk" "$INSTDIR\uninst.exe"
SectionEnd

Section -Post
  WriteUninstaller "$INSTDIR\uninst.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayName" "Percolator"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString" "$INSTDIR\uninst.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayVersion" "${PRODUCT_VERSION}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "URLInfoAbout" "${PRODUCT_WEB_SITE}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "Publisher" "${PRODUCT_PUBLISHER}"
  WriteRegStr HKLM "Software/Microsoft/Windows/CurrentVersion/Uninstall/Product" "${PRODUCT_NAME} ${PRODUCT_VERSION}" "$INSTDIR\uninst.exe"
  
  ;Exec "$INSTDIR\dbUpgrader.exe"  ; this would run a post-install program you provided
SectionEnd

Function un.onUninstSuccess
  HideWindow
  MessageBox MB_ICONINFORMATION|MB_OK "Percolator was successfully removed from your computer."
FunctionEnd

Function un.onInit
  MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove Percolator and all of its components?" IDYES +2
  Abort
FunctionEnd

Function .onInit
  SetOutPath $TEMP
  File /oname=setup.bmp "${SETUP_BITMAP}"

  advsplash::show 1500 500 500 -1 $TEMP\setup

  Pop $0 ; $0 has '1' if the user closed the splash screen early,
         ; '0' if everything closed normal, and '-1' if some error occured.

  Delete $TEMP\setup.bmp
FunctionEnd

Section Uninstall
;  UnRegDLL "yourRegistered.dll"

  Delete "$INSTDIR\${PRODUCT_NAME}.url"
  Delete "$INSTDIR\uninst.exe"
  Delete "$INSTDIR\readme.txt"
  Delete "${MING_BIN}\boost_filesystem-gcc45-mt-1_46_1.dll"
  Delete "${MING_BIN}\boost_system-gcc45-mt-1_46_1.dll"
  Delete "${MING_BIN}\libgcc_s_sjlj-1.dll"
  Delete "${MING_BIN}\libstdc++-6.dll"
  Delete "${MING_BIN}\libxerces-c-3-0.dll"
  Delete "$INSTDIR\${PROGICON}"
  Delete "$INSTDIR\license.txt"
;  Delete "$INSTDIR\dbUpgrader.exe"

  SetShellVarContext all
  Delete "$SMPROGRAMS\${PRODUCT_NAME}\${PRODUCT_NAME}.lnk"
  Delete "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall.lnk"
  Delete "$SMPROGRAMS\${PRODUCT_NAME}\Website.lnk"
  Delete "$SMPROGRAMS\${PRODUCT_NAME}\Help.lnk"
  
  Delete "$DESKTOP\${PRODUCT_NAME}.lnk"

  RMDir "$SMPROGRAMS\${PRODUCT_NAME}"
  RMDir "$INSTDIR"

  DeleteRegKey ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}"
  SetAutoClose true
SectionEnd