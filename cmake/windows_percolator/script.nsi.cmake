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

;-----------------------------------------------------------------------------
; Installer version
;-----------------------------------------------------------------------------

!define VER_MAJOR "@CPACK_PACKAGE_VERSION_MAJOR@"
!define VER_MINOR "@CPACK_PACKAGE_VERSION_MINOR@"
!define VER_BUILD "@CPACK_PACKAGE_VERSION_PATCH@"
!define VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"
!define REVISION ""


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
!define MUI_WELCOMEFINISHPAGE_BITMAP ${NSI_PATH}\welcome.bmp
!define MUI_WELCOMEPAGE_TITLE "@CPACK_PACKAGE_NAME@ ${VERSION} Setup$\r$\nInstaller"
!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation.$\r$\n$\r$\n$_CLICK"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP ${NSI_PATH}\page_header.bmp
!define MUI_COMPONENTSPAGE_SMALLDESC
!define MUI_FINISHPAGE_TITLE "Percolator Completed"
!define MUI_FINISHPAGE_LINK "Click here to visit the @CPACK_PACKAGE_NAME@ website."
!define MUI_FINISHPAGE_LINK_LOCATION "http://per-colator.com/"
!define MUI_FINISHPAGE_NOAUTOCLOSE
!define MUI_FINISHPAGE_RUN
!define MUI_FINISHPAGE_RUN_FUNCTION "LaunchPercolator"

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


Name "percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@"
Caption "Percolator Installer"
outfile "@CMAKE_BINARY_DIR@/percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@-win32.exe"
installDir $PROGRAMFILES\Percolator
InstallDirRegKey HKCU "Software\Percolator" ""
ShowInstDetails show
ShowUnInstDetails show
ReserveFile NSIS.InstallOptions.ini
ReserveFile "${NSISDIR}\Plugins\InstallOptions.dll"

##############################################################################
# #
# FINISH PAGE LAUNCHER FUNCTIONS #
# #
##############################################################################

Function LaunchPercolator
   Exec "$INSTDIR\percolator.exe"
FunctionEnd

##############################################################################
# #
# RE-INSTALLER FUNCTIONS #
# #
##############################################################################

Function PageReinstall
   ReadRegStr $R0 HKLM "Software\Percolator" ""
   StrCmp $R0 "" 0 +2
   Abort

   ;Detect version
   ReadRegDWORD $R0 HKLM "Software\Percolator" "VersionMajor"
   IntCmp $R0 ${VER_MAJOR} minor_check new_version older_version
   minor_check:
      ReadRegDWORD $R0 HKLM "Software\Percolator" "VersionMinor"
      IntCmp $R0 ${VER_MINOR} build_check new_version older_version
   build_check:
      ReadRegDWORD $R0 HKLM "Software\Percolator" "VersionBuild"
      IntCmp $R0 ${VER_BUILD} same_version new_version older_version

   new_version:
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 1" "Text" "An older version of Percolator is installed on your system. It is recommended that you uninstall the current version before installing. Select the operation you want to perform and click Next to continue."
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 2" "Text" "Uninstall before installing"
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 3" "Text" "Do not uninstall"
      !insertmacro MUI_HEADER_TEXT "Already Installed" "Choose how you want to install Percolator."
      StrCpy $R0 "1"
      Goto reinst_start

   older_version:
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 1" "Text" "A newer version of Percolator is already installed! It is not recommended that you install an older version. If you really want to install this older version, it is better to uninstall the current version first. Select the operation you want to perform and click Next to continue."
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 2" "Text" "Uninstall before installing"
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 3" "Text" "Do not uninstall"
      !insertmacro MUI_HEADER_TEXT "Already Installed" "Choose how you want to install Percolator."
      StrCpy $R0 "1"
      Goto reinst_start

   same_version:
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 1" "Text" "Percolator ${VERSION} is already installed.\r\nSelect the operation you want to perform and click Next to continue."
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 2" "Text" "Add/Reinstall components"
      !insertmacro INSTALLOPTIONS_WRITE "NSIS.InstallOptions.ini" "Field 3" "Text" "Uninstall Percolator"
      !insertmacro MUI_HEADER_TEXT "Already Installed" "Choose the maintenance option to perform."
      StrCpy $R0 "2"

   reinst_start:
      !insertmacro INSTALLOPTIONS_DISPLAY "NSIS.InstallOptions.ini"
FunctionEnd

Function PageLeaveReinstall
   !insertmacro INSTALLOPTIONS_READ $R1 "NSIS.InstallOptions.ini" "Field 2" "State"
   StrCmp $R0 "1" 0 +2
   StrCmp $R1 "1" reinst_uninstall reinst_done
   StrCmp $R0 "2" 0 +3
   StrCmp $R1 "1" reinst_done reinst_uninstall
   reinst_uninstall:
      ReadRegStr $R1 HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "UninstallString"
      HideWindow
      ClearErrors
      ExecWait '$R1 _?=$INSTDIR'
      IfErrors no_remove_uninstaller
      IfFileExists "$INSTDIR\percolator.exe" no_remove_uninstaller
      Delete $R1
      RMDir $INSTDIR
   no_remove_uninstaller:
      StrCmp $R0 "2" 0 +3
      UAC::Unload
      Quit
      BringToFront
   reinst_done:
FunctionEnd

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
SectionEnd

Section "Pin Xml Schema" SecPinXml
  SetOutPath "@PIN_SCHEMA_LOCATION@"
  file @CMAKE_SOURCE_DIR@/src/xml/percolator_in.xsd
SectionEnd

Section "Pout Xml Schema" SecPoutXml
  SetOutPath "@POUT_SCHEMA_LOCATION@"
  file @CMAKE_SOURCE_DIR@/src/xml/percolator_out.xsd
SectionEnd

Section "Start Menu Program Group" SEC_START_MENU
  SectionIn 1 2
  SetDetailsPrint textonly
  DetailPrint "Adding shortcuts for the Percolator program group to the Start Menu."
  SetDetailsPrint listonly
  SetShellVarContext all
  RMDir /r "$SMPROGRAMS\Percolator"
  CreateDirectory "$SMPROGRAMS\Percolator"
# CreateShortCut "$SMPROGRAMS\Percolator\LICENSE.lnk" "$INSTDIR\LICENSE.txt"
  CreateShortCut "$SMPROGRAMS\Percolator\Percolator.lnk" "$INSTDIR\percolator.exe"
# CreateShortCut "$SMPROGRAMS\Percolator\Release notes.lnk" "$INSTDIR\NOTES.txt"
  CreateShortCut "$SMPROGRAMS\Percolator\Uninstall.lnk" "$INSTDIR\uninstall.exe"
  SetShellVarContext current
SectionEnd

Section "Desktop Shortcut" SEC_DESKTOP
      SectionIn 1 2
      SetDetailsPrint textonly
      DetailPrint "Creating Desktop Shortcuts"
      SetDetailsPrint listonly
      CreateShortCut "$DESKTOP\Percolator.lnk" "$INSTDIR\percolator.exe"
SectionEnd

Section -AdditionalIcons
  WriteIniStr "$INSTDIR\${PRODUCT_NAME}.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Website.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  SetShellVarContext all		; scope is "All Users"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall.lnk" "$INSTDIR\uninst.exe"
SectionEnd


; Installer section descriptions
;--------------------------------
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
!insertmacro MUI_DESCRIPTION_TEXT ${SEC01} "Percolator essentials."
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_START_MENU} "Percolator program group."
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_DESKTOP} "Desktop shortcut for Percolator."
!insertmacro MUI_FUNCTION_DESCRIPTION_END

Section -Post
   ;Uninstaller file.
   SetDetailsPrint textonly
   DetailPrint "Writing Uninstaller"
   SetDetailsPrint listonly
   WriteUninstaller $INSTDIR\uninstall.exe

   ;Registry keys required for installer version handling and uninstaller.
   SetDetailsPrint textonly
   DetailPrint "Writing Installer Registry Keys"
   SetDetailsPrint listonly

   ;Version numbers used to detect existing installation version for comparisson.
   WriteRegStr HKLM "Software\Percolator" "" $INSTDIR
   WriteRegDWORD HKLM "Software\Percolator" "VersionMajor" "${VER_MAJOR}"
   WriteRegDWORD HKLM "Software\Percolator" "VersionMinor" "${VER_MINOR}"
   WriteRegDWORD HKLM "Software\Percolator" "VersionRevision" "${REVISION}"
   WriteRegDWORD HKLM "Software\Percolator" "VersionBuild" "${VER_BUILD}"

   ;Add or Remove Programs entry.
   WriteRegExpandStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "UninstallString" "$INSTDIR\Uninstall.exe"
   WriteRegExpandStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "InstallLocation" "$INSTDIR"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "DisplayName" "Percolator"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "Publisher" "http://per-colator.com/"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "DisplayIcon" "$INSTDIR\Uninstall.exe,0"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "DisplayVersion" "${VERSION}"
   WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "VersionMajor" "${VER_MAJOR}"
   WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "VersionMinor" "${VER_MINOR}.${REVISION}"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "URLInfoAbout" "http://per-colator.com/"
   WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "HelpLink" "http://per-colator.com/"
   WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "NoModify" "1"
   WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "NoRepair" "1"

   ; Register Percolator:// protocol handler
   WriteRegStr HKCR "percolator" "" "URL: Percolator Protocol"
   WriteRegStr HKCR "percolator\DefaultIcon" "" $INSTDIR\percolator.exe,1
   WriteRegStr HKCR "percolator\shell" "" "open"
   WriteRegStr HKCR "percolator\shell\open\command" "" "$INSTDIR\percolator.exe %1"

   SetDetailsPrint textonly
   DetailPrint "Finsihed."
SectionEnd

##############################################################################
# #
# UNINSTALLER SECTION #
# #
##############################################################################

Var UnPageUserAppDataDialog
Var UnPageUserAppDataCheckbox
Var UnPageUserAppDataCheckbox_State
Var UnPageUserAppDataEditBox

Function un.UnPageUserAppData
   !insertmacro MUI_HEADER_TEXT "Uninstall Percolator" "Remove Percolator's data folder from your computer."
   nsDialogs::Create /NOUNLOAD 1018
   Pop $UnPageUserAppDataDialog

   ${If} $UnPageUserAppDataDialog == error
      Abort
   ${EndIf}

   ${NSD_CreateLabel} 0 0 100% 12u "Do you want to delete Percolator's data folder?"
   Pop $0

   ${NSD_CreateText} 0 13u 100% 12u "$LOCALAPPDATA\Percolator"
   Pop $UnPageUserAppDataEditBox
   SendMessage $UnPageUserAppDataEditBox ${EM_SETREADONLY} 1 0

   ${NSD_CreateLabel} 0 46u 100% 24u "Leave unchecked to keep the data folder for later use or check to delete the data folder."
   Pop $0

   ${NSD_CreateCheckbox} 0 71u 100% 8u "Yes, delete this data folder."
   Pop $UnPageUserAppDataCheckbox

   nsDialogs::Show
FunctionEnd

Function un.UnPageUserAppDataLeave
   ${NSD_GetState} $UnPageUserAppDataCheckbox $UnPageUserAppDataCheckbox_State
FunctionEnd

Section Uninstall
   IfFileExists "$INSTDIR\percolator.exe" percolator_installed
      MessageBox MB_YESNO "It does not appear that Percolator is installed in the directory '$INSTDIR'.$\r$\nContinue anyway (not recommended)?" IDYES percolator_installed
      Abort "Uninstall aborted by user"
   percolator_installed:

   ;Delete registry keys.
   DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator"
   DeleteRegValue HKLM "Software\Percolator" "VersionBuild"
   DeleteRegValue HKLM "Software\Percolator" "VersionMajor"
   DeleteRegValue HKLM "Software\Percolator" "VersionMinor"
   DeleteRegValue HKLM "Software\Percolator" "VersionRevision"
   DeleteRegValue HKLM "Software\Percolator" ""
   DeleteRegKey HKLM "Software\Percolator"

   DeleteRegKey HKCR "percolator"

   ;Start menu shortcuts.
   !ifdef OPTION_SECTION_SC_START_MENU
      SetShellVarContext all
      RMDir /r "$SMPROGRAMS\Percolator"
      SetShellVarContext current
   !endif

   ;Desktop shortcut.
   !ifdef OPTION_SECTION_SC_DESKTOP
      IfFileExists "$DESKTOP\Percolator.lnk" 0 +2
         Delete "$DESKTOP\Percolator.lnk"
   !endif

   ;Quick Launch shortcut.
   !ifdef OPTION_SECTION_SC_QUICK_LAUNCH
      IfFileExists "$QUICKLAUNCH\Percolator.lnk" 0 +2
         Delete "$QUICKLAUNCH\Percolator.lnk"
   !endif

   ;Remove all the Program Files.
   RMDir /r $INSTDIR

   ;Uninstall User Data if option is checked, otherwise skip.
   ${If} $UnPageUserAppDataCheckbox_State == ${BST_CHECKED}
      RMDir /r "$LOCALAPPDATA\Percolator"
   ${EndIf}

   SetDetailsPrint textonly
   DetailPrint "Finsihed."
SectionEnd

##############################################################################
# #
# NSIS Installer Event Handler Functions #
# #
##############################################################################

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
  Delete "$INSTDIR\boost_filesystem-gcc45-mt-1_46_1.dll"
  Delete "$INSTDIR\boost_system-gcc45-mt-1_46_1.dll"
  Delete "$INSTDIR\libgcc_s_sjlj-1.dll"
  Delete "$INSTDIR\libstdc++-6.dll"
  Delete "$INSTDIR\libxerces-c-3-0.dll"
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