; Percolator installer

;-----------------------------------------------------------------------------
; Some paths.
;-----------------------------------------------------------------------------

!define MING_BIN "@MING_PATH@/bin"
!define MING_LIB "@MING_PATH@/lib"
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

;----------------------------------------------------------------------------------------------
; Installer configuration
;-------------------------------------------------------------------------------------------

; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "Percolator"
!define PRODUCT_VERSION "Percolator-@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"
!define PRODUCT_PUBLISHER "<Lukas Kall>"
!define PRODUCT_WEB_SITE "http://per-colator.com/"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
!define PRODUCT_UNINST_ROOT_KEY "HKLM"
!define EXECUTABLE "percolator.exe"
!define PROGICON "${NSI_PATH}\installer.ico"
!define PROGICON2 "${NSI_PATH}\uninstall.ico"
!define PROGICON3 "${NSI_PATH}\percolator.ico"
!define SETUP_BITMAP "${NSI_PATH}\welcome.bmp"


; MUI 1.67 compatible ------
!include "MUI2.nsh"
!include LogicLib.nsh ;Used by APPDATA uninstaller.
!include nsDialogs.nsh ;Used by APPDATA uninstaller.
!include InstallOptions.nsh

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_LANGDLL_ALLLANGUAGES
;--------------------------------
;Language Selection Dialog Settings
;Remember the installer language
!define MUI_LANGDLL_REGISTRY_ROOT "HKCU" 
!define MUI_LANGDLL_REGISTRY_KEY "Software\Modern UI Test" 
!define MUI_LANGDLL_REGISTRY_VALUENAME "Installer Language"
!define MUI_ICON ${PROGICON}
!define MUI_UNICON ${PROGICON2}  ; define uninstall icon to appear in Add/Remove Programs
!define MUI_WELCOMEFINISHPAGE_BITMAP ${NSI_PATH}\welcome.bmp
!define MUI_WELCOMEPAGE_TITLE "Percolator - @CPACK_PACKAGE_NAME@ ${VERSION} Setup$\r$\nInstaller"
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
Page custom PageReinstall PageLeaveReinstall
# ; License page
# !insertmacro MUI_PAGE_LICENSE "license.txt"   ; text file with license terms
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\readme.txt"  ; readme.txt file for user
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_INSTFILES
;--------------------------------
;Languages

  !insertmacro MUI_LANGUAGE "English" ;first language is the default language
  !insertmacro MUI_LANGUAGE "French"
  !insertmacro MUI_LANGUAGE "German"
  !insertmacro MUI_LANGUAGE "Spanish"
  !insertmacro MUI_LANGUAGE "SpanishInternational"
  !insertmacro MUI_LANGUAGE "SimpChinese"
  !insertmacro MUI_LANGUAGE "TradChinese"
  !insertmacro MUI_LANGUAGE "Japanese"
  !insertmacro MUI_LANGUAGE "Korean"
  !insertmacro MUI_LANGUAGE "Italian"
  !insertmacro MUI_LANGUAGE "Dutch"
  !insertmacro MUI_LANGUAGE "Danish"
  !insertmacro MUI_LANGUAGE "Swedish"
  !insertmacro MUI_LANGUAGE "Norwegian"
  !insertmacro MUI_LANGUAGE "NorwegianNynorsk"
  !insertmacro MUI_LANGUAGE "Finnish"
  !insertmacro MUI_LANGUAGE "Greek"
  !insertmacro MUI_LANGUAGE "Russian"
  !insertmacro MUI_LANGUAGE "Portuguese"
  !insertmacro MUI_LANGUAGE "PortugueseBR"
  !insertmacro MUI_LANGUAGE "Polish"
  !insertmacro MUI_LANGUAGE "Ukrainian"
  !insertmacro MUI_LANGUAGE "Czech"
  !insertmacro MUI_LANGUAGE "Slovak"
  !insertmacro MUI_LANGUAGE "Croatian"
  !insertmacro MUI_LANGUAGE "Bulgarian"
  !insertmacro MUI_LANGUAGE "Hungarian"
  !insertmacro MUI_LANGUAGE "Thai"
  !insertmacro MUI_LANGUAGE "Romanian"
  !insertmacro MUI_LANGUAGE "Latvian"
  !insertmacro MUI_LANGUAGE "Macedonian"
  !insertmacro MUI_LANGUAGE "Estonian"
  !insertmacro MUI_LANGUAGE "Turkish"
  !insertmacro MUI_LANGUAGE "Lithuanian"
  !insertmacro MUI_LANGUAGE "Slovenian"
  !insertmacro MUI_LANGUAGE "Serbian"
  !insertmacro MUI_LANGUAGE "SerbianLatin"
  !insertmacro MUI_LANGUAGE "Arabic"
  !insertmacro MUI_LANGUAGE "Farsi"
  !insertmacro MUI_LANGUAGE "Hebrew"
  !insertmacro MUI_LANGUAGE "Indonesian"
  !insertmacro MUI_LANGUAGE "Mongolian"
  !insertmacro MUI_LANGUAGE "Luxembourgish"
  !insertmacro MUI_LANGUAGE "Albanian"
  !insertmacro MUI_LANGUAGE "Breton"
  !insertmacro MUI_LANGUAGE "Belarusian"
  !insertmacro MUI_LANGUAGE "Icelandic"
  !insertmacro MUI_LANGUAGE "Malay"
  !insertmacro MUI_LANGUAGE "Bosnian"
  !insertmacro MUI_LANGUAGE "Kurdish"
  !insertmacro MUI_LANGUAGE "Irish"
  !insertmacro MUI_LANGUAGE "Uzbek"
  !insertmacro MUI_LANGUAGE "Galician"
  !insertmacro MUI_LANGUAGE "Afrikaans"
  !insertmacro MUI_LANGUAGE "Catalan"
  !insertmacro MUI_LANGUAGE "Esperanto"


;--------------------------------
;Reserve Files
  
  ;If you are using solid compression, files that are required before
  ;the actual installation should be stored first in the data block,
  ;because this will make your installer start faster.
  
  !insertmacro MUI_RESERVEFILE_LANGDLL

;--------------------------------
UninstPage custom un.UnPageUserAppData un.UnPageUserAppDataLeave
# !insertmacro MUI_LANGUAGE "English"

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
   File "${BUILD_PATH}\src\qvality\qvality.exe"
   File "${NSI_PATH}\percolator.ico"
   ;MINGW
   File "${BUILD_PATH}\\nsi_bin\\*.*"

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
  CreateShortCut "$SMPROGRAMS\Percolator\Percolator.lnk" "$INSTDIR\qvality.exe"
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
      CreateShortCut "$DESKTOP\Percolator.lnk" "$INSTDIR\qvality.exe"
SectionEnd

Section "Additional Icons" SEC_ADD_ICONS
  WriteIniStr "$INSTDIR\${PRODUCT_NAME}.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Website.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  SetShellVarContext all		; scope is "All Users"
  CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall.lnk" "$INSTDIR\uninst.exe"
SectionEnd

# LangString Desc_Section1 ${LANG_ENGLISH} "Percolator essentials."
# LangString Desc_Section2 ${LANG_ENGLISH} "Percolator program group."
# LangString Desc_Section3 ${LANG_ENGLISH} "Desktop shortcut for Percolator."
# LangString Desc_Section4 ${LANG_ENGLISH} "Additional Icons for Percolator."
# 
# ; Installer section descriptions
# ;--------------------------------
# !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
#   !insertmacro MUI_DESCRIPTION_TEXT ${SEC01} ${Desc_Section1}
#   !insertmacro MUI_DESCRIPTION_TEXT ${SEC_START_MENU} ${Desc_Section2}
#   !insertmacro MUI_DESCRIPTION_TEXT ${SEC_DESKTOP} ${Desc_Section3}
#   !insertmacro MUI_DESCRIPTION_TEXT ${SEC_ADD_ICONS} ${Desc_Section4}
# !insertmacro MUI_FUNCTION_DESCRIPTION_END

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
   WriteRegStr HKCR "percolator\DefaultIcon" "" $INSTDIR\qvality.exe,1
   WriteRegStr HKCR "percolator\shell" "" "open"
   WriteRegStr HKCR "percolator\shell\open\command" "" "$INSTDIR\percolator.exe %1"
   WriteRegStr HKCR "percolator\shell\open\command" "" "$INSTDIR\qvality.exe %1"

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
      IfFileExists "$DESKTOP\qvality.lnk" 0 +2
         Delete "$DESKTOP\qvality.lnk"
   !endif

   ;Quick Launch shortcut.
   !ifdef OPTION_SECTION_SC_QUICK_LAUNCH
      IfFileExists "$QUICKLAUNCH\Percolator.lnk" 0 +2
         Delete "$QUICKLAUNCH\Percolator.lnk"
      IfFileExists "$QUICKLAUNCH\qvality.lnk" 0 +2
         Delete "$QUICKLAUNCH\qvality.lnk"
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
  !insertmacro MUI_UNGETLANGUAGE
  MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove Percolator and all of its components?" IDYES +2
  Abort
FunctionEnd

;Function that calls a messagebox when installation finished correctly
Function .onInstSuccess
  MessageBox MB_OK "You have successfully installed Percolator. Use the desktop icon to start the program."
FunctionEnd

Function .onInit
  !insertmacro MUI_LANGDLL_DISPLAY
  SetOutPath $TEMP
  File /oname=setup.bmp "${SETUP_BITMAP}"

  advsplash::show 1500 500 500 -1 $TEMP\setup

  Pop $0 ; $0 has '1' if the user closed the splash screen early,
         ; '0' if everything closed normal, and '-1' if some error occured.

  Delete $TEMP\setup.bmp
FunctionEnd

Section Uninstall
  Delete "$INSTDIR\${PRODUCT_NAME}.url"
  Delete "$INSTDIR\uninst.exe"
  Delete "$INSTDIR\readme.txt"
  Delete "$INSTDIR\boost_filesystem-gcc45-mt-1_46_1.dll"
  Delete "$INSTDIR\boost_system-gcc45-mt-1_46_1.dll"
  Delete "$INSTDIR\libgcc_s_sjlj-1.dll"
  Delete "$INSTDIR\libstdc++-6.dll"
  Delete "$INSTDIR\libxerces-c-3-0.dll"
  Delete "$INSTDIR\${PROGICON}"
  Delete "$INSTDIR\${PROGICON2}"
  Delete "$INSTDIR\${PROGICON3}"
  Delete "$INSTDIR\license.txt"

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
