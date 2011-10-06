!include "MUI2.nsh"

;Percolator installer script.

;-----------------------------------------------------------------------------
; Some installer script options (comment-out options not required)
;-----------------------------------------------------------------------------
;!define OPTION_LICENSE_AGREEMENT
!define OPTION_UAC_PLUGIN_ENHANCED
!define OPTION_SECTION_SC_START_MENU
!define OPTION_SECTION_SC_DESKTOP
!define OPTION_SECTION_SC_QUICK_LAUNCH
!define OPTION_FINISHPAGE
!define OPTION_FINISHPAGE_LAUNCHER
!define OPTION_FINISHPAGE_RELEASE_NOTES


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
!addplugindir ${NSI_PATH}
!addplugindir ${NSI_PATH}/nsis_uac
!addplugindir ${NSI_PATH}/nsis_processes
!addplugindir ${NSI_PATH}/nsis_processes/src
;-----------------------------------------------------------------------------
; Installer version
;-----------------------------------------------------------------------------

!define VER_MAJOR "@CPACK_PACKAGE_VERSION_MAJOR@"
!define VER_MINOR "@CPACK_PACKAGE_VERSION_MINOR@"
!define VER_BUILD "@CPACK_PACKAGE_VERSION_PATCH@"
!define VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"

;-----------------------------------------------------------------------------
; Installer build timestamp.
;-----------------------------------------------------------------------------
!define /date BUILD_TIME "built on %Y/%m/%d at %I:%M %p"

;-----------------------------------------------------------------------------
; Initial installer setup and definitions.
;-----------------------------------------------------------------------------
Name "percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@"
Caption "Percolator Installer"
outfile "@CMAKE_BINARY_DIR@/percolator-@CPACK_PACKAGE_VERSION_MAJOR@_@CPACK_PACKAGE_VERSION_MINOR@-win32.exe"
installDir $PROGRAMFILES\Percolator
InstallDirRegKey HKCU "Software\Percolator" ""
InstType Standard
InstType Full
InstType Minimal
CRCCheck On
RequestExecutionLevel user ;Now using the UAC plugin.
ReserveFile NSIS.InstallOptions.ini
ReserveFile "${NSISDIR}\Plugins\InstallOptions.dll"


@CPACK_NSIS_SECTION_SELECTED_VARS@
;-----------------------------------------------------------------------------
; Include some required header files.
;-----------------------------------------------------------------------------
!include LogicLib.nsh ;Used by APPDATA uninstaller.
!include nsDialogs.nsh ;Used by APPDATA uninstaller.
!include MUI2.nsh ;Used by APPDATA uninstaller.
!include InstallOptions.nsh ;Required by MUI2 to support old MUI_INSTALLOPTIONS.
!include Memento.nsh ;Remember user selections.
!include WinVer.nsh ;Windows version detection.
!include WordFunc.nsh ;Used by VersionCompare macro function.
!include ${NSI_PATH}\nsis_uac\UAC.nsh ;Used by the UAC elevation to install as user or admin.
!include ${NSI_PATH}/nsis_processes/src/ProcFunc.nsh 
;-----------------------------------------------------------------------------
; Memento selections stored in registry.
;-----------------------------------------------------------------------------
!define MEMENTO_REGISTRY_ROOT HKLM
!define MEMENTO_REGISTRY_KEY Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator


;-----------------------------------------------------------------------------
; Modern User Interface (MUI) defintions and setup.
;-----------------------------------------------------------------------------

!define MUI_ABORTWARNING
!define MUI_ICON ${NSI_PATH}\installer.ico
!define MUI_UNICON ${NSI_PATH}\installer.ico
!define MUI_WELCOMEFINISHPAGE_BITMAP ${NSI_PATH}\welcome.bmp
!define MUI_WELCOMEPAGE_TITLE "@CPACK_PACKAGE_NAME@ ${VERSION} Setup$\r$\nInstaller"
!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation.$\r$\n$\r$\n$_CLICK"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP ${NSI_PATH}\page_header.bmp
!define MUI_COMPONENTSPAGE_SMALLDESC
!define MUI_FINISHPAGE_TITLE "@CPACK_PACKAGE_NAME@ Install Completed"
!define MUI_FINISHPAGE_LINK "Click here to visit the @CPACK_PACKAGE_NAME@ website."
!define MUI_FINISHPAGE_LINK_LOCATION "http://per-colator.com/"
!define MUI_FINISHPAGE_NOREBOOTSUPPORT
!ifdef OPTION_FINISHPAGE_RELEASE_NOTES
   !define MUI_FINISHPAGE_SHOWREADME_NOTCHECKED
   !define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\NOTES.txt"
   !define MUI_FINISHPAGE_SHOWREADME_TEXT "Show release notes"
!endif
!ifdef OPTION_FINISHPAGE_LAUNCHER
   !define MUI_FINISHPAGE_NOAUTOCLOSE
   !define MUI_FINISHPAGE_RUN
   !define MUI_FINISHPAGE_RUN_FUNCTION "LaunchPercolator"
!endif


;-----------------------------------------------------------------------------
; Page macros.
;-----------------------------------------------------------------------------

!insertmacro MUI_PAGE_LICENSE "@CMAKE_SOURCE_DIR@/src/COPYING"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Page custom PageReinstall PageLeaveReinstall
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!ifdef OPTION_FINISHPAGE
   !insertmacro MUI_PAGE_FINISH
!endif
!insertmacro MUI_UNPAGE_CONFIRM
UninstPage custom un.UnPageUserAppData un.UnPageUserAppDataLeave
!insertmacro MUI_UNPAGE_INSTFILES

##############################################################################
# #
# PROCESS HANDLING FUNCTIONS AND MACROS #
# #
##############################################################################

!macro CheckForProcess processName gotoWhenFound gotoWhenNotFound
   FindProcess ${processName}
   StrCmp $R0 "0" ${gotoWhenNotFound} ${gotoWhenFound}
!macroend

!macro ConfirmEndProcess processName
   MessageBox MB_YESNO|MB_ICONEXCLAMATION \
     "Found ${processName} process(s) which need to be stopped.$\nDo you want the installer to stop these for you?" \
     IDYES process_${processName}_kill IDNO process_${processName}_ended
   process_${processName}_kill:
      DetailPrint "Killing ${processName} processes."
      KillProcess ${processName}
      Sleep 1500
      StrCmp $R0 "1" process_${processName}_ended
      DetailPrint "Process to kill not found!"
   process_${processName}_ended:
!macroend

!macro CheckAndConfirmEndProcess processName
   !insertmacro CheckForProcess ${processName} 0 no_process_${processName}_to_end
   !insertmacro ConfirmEndProcess ${processName}
   no_process_${processName}_to_end:
!macroend

Function EnsurePercolatorShutdown
   !insertmacro CheckAndConfirmEndProcess "percolator.exe"
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

##############################################################################
# #
# INSTALLER SECTIONS #
# #
##############################################################################

Section "Dummy Section" SecDummy
   SectionIn 1 2 3 RO
   SetDetailsPrint listonly

   SetDetailsPrint textonly
   DetailPrint "Installing Percolator essentials."
   SetDetailsPrint listonly
   SetOutPath "$INSTDIR"
   
   !ifdef INSTALL_PATH
        ;Main executable.
        File "${INSTALL_PATH}\bin\percolator.exe"
   !endif
   !ifndef INSTALL_PATH
        ;Main executable.
        File "${BUILD_PATH}\percolator.exe"
        File @CMAKE_BINARY_DIR@/src/percolator.exe
   !endif
   
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

SectionGroup "Shortcuts"

!ifdef OPTION_SECTION_SC_START_MENU
   ${MementoSection} "Start Menu Program Group" SEC_START_MENU
      SectionIn 1 2
      SetDetailsPrint textonly
      DetailPrint "Adding shortcuts for the Percolator program group to the Start Menu."
      SetDetailsPrint listonly
      SetShellVarContext all
      RMDir /r "$SMPROGRAMS\Percolator"
      CreateDirectory "$SMPROGRAMS\Percolator"
      CreateShortCut "$SMPROGRAMS\Percolator\LICENSE.lnk" "$INSTDIR\LICENSE.txt"
      CreateShortCut "$SMPROGRAMS\Percolator\Percolator.lnk" "$INSTDIR\percolator.exe"
      CreateShortCut "$SMPROGRAMS\Percolator\Release notes.lnk" "$INSTDIR\NOTES.txt"
      CreateShortCut "$SMPROGRAMS\Percolator\Uninstall.lnk" "$INSTDIR\uninstall.exe"
      SetShellVarContext current
   ${MementoSectionEnd}
!endif

!ifdef OPTION_SECTION_SC_DESKTOP
   ${MementoSection} "Desktop Shortcut" SEC_DESKTOP
      SectionIn 1 2
      SetDetailsPrint textonly
      DetailPrint "Creating Desktop Shortcuts"
      SetDetailsPrint listonly
      CreateShortCut "$DESKTOP\Percolator.lnk" "$INSTDIR\percolator.exe"
   ${MementoSectionEnd}
!endif

!ifdef OPTION_SECTION_SC_QUICK_LAUNCH
   ${MementoSection} "Quick Launch Shortcut" SEC_QUICK_LAUNCH
      SectionIn 1 2
      SetDetailsPrint textonly
      DetailPrint "Creating Quick Launch Shortcut"
      SetDetailsPrint listonly
      CreateShortCut "$QUICKLAUNCH\Percolator.lnk" "$INSTDIR\percolator.exe"
   ${MementoSectionEnd}
!endif

SectionGroupEnd


; Installer section descriptions
;--------------------------------
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_PERCOLATOR} "Percolator essentials."
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_START_MENU} "Percolator program group."
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_DESKTOP} "Desktop shortcut for Percolator."
!insertmacro MUI_DESCRIPTION_TEXT ${SEC_QUICK_LAUNCH} "Quick Launch shortcut for Percolator."
!insertmacro MUI_FUNCTION_DESCRIPTION_END

Section -post

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
   WriteRegExpandStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "UninstallString" '"$INSTDIR\Uninstall.exe"'
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
   WriteRegStr HKCR "percolator\shell\open\command" "" '"$INSTDIR\percolator.exe" "%1"'

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

Function .onInit
   !insertmacro INSTALLOPTIONS_EXTRACT "NSIS.InstallOptions.ini"

   ;Remove Quick Launch option from Windows 7, as no longer applicable - usually.
   ${IfNot} ${AtMostWinVista}
      SectionSetText ${SEC_QUICK_LAUNCH} "Quick Launch Shortcut (N/A)"
      SectionSetFlags ${SEC_QUICK_LAUNCH} ${SF_RO}
      SectionSetInstTypes ${SEC_QUICK_LAUNCH} 0
   ${EndIf}

   ${MementoSectionRestore}

   UAC_Elevate:
      UAC::RunElevated
      StrCmp 1223 $0 UAC_ElevationAborted ; UAC dialog aborted by user?
      StrCmp 0 $0 0 UAC_Err ; Error?
      StrCmp 1 $1 0 UAC_Success ;Are we the real deal or just the wrapper?
      Quit

   UAC_Err:
      MessageBox MB_ICONSTOP "Unable to elevate, error $0"
      Abort

   UAC_ElevationAborted:
      Abort

   UAC_Success:
      StrCmp 1 $3 +4 ;Admin?
      StrCmp 3 $1 0 UAC_ElevationAborted ;Try again?
      MessageBox MB_ICONSTOP "This installer requires admin access, try again"
      goto UAC_Elevate

   ;Prevent multiple instances.
   System::Call 'kernel32::CreateMutexA(i 0, i 0, t "percolatorInstaller") i .r1 ?e'
   Pop $R0
   StrCmp $R0 0 +3
      MessageBox MB_OK|MB_ICONEXCLAMATION "The installer is already running."
      Abort

   ;Use available InstallLocation when possible. This is useful in the uninstaller
   ;via re-install, which would otherwise use a default location - a bug.
   ReadRegStr $R0 HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Percolator" "InstallLocation"
   StrCmp $R0 "" SkipSetInstDir
   StrCpy $INSTDIR $R0
   SkipSetInstDir:

   ;Shutdown percolator in case Add/Remove re-installer option used.
   Call EnsurePercolatorkShutdown
FunctionEnd

Function .onInstSuccess
   ${MementoSectionSave}
   UAC::Unload ;Must call unload!
FunctionEnd

Function .onInstFailed
   UAC::Unload ;Must call unload!
FunctionEnd

##############################################################################
# #
# NSIS Uninstaller Event Handler Functions #
# #
##############################################################################

Function un.onInit
   UAC_Elevate:
      UAC::RunElevated
      StrCmp 1223 $0 UAC_ElevationAborted ; UAC dialog aborted by user?
      StrCmp 0 $0 0 UAC_Err ; Error?
      StrCmp 1 $1 0 UAC_Success ;Are we the real deal or just the wrapper?
      Quit

   UAC_Err:
      MessageBox MB_ICONSTOP "Unable to elevate, error $0"
      Abort

   UAC_ElevationAborted:
      Abort

   UAC_Success:
      StrCmp 1 $3 +4 ;Admin?
      StrCmp 3 $1 0 UAC_ElevationAborted ;Try again?
      MessageBox MB_ICONSTOP "This uninstaller requires admin access, try again"
      goto UAC_Elevate

   ;Prevent multiple instances.
   System::Call 'kernel32::CreateMutexA(i 0, i 0, t "percolatorUninstaller") i .r1 ?e'
   Pop $R0
   StrCmp $R0 0 +3
      MessageBox MB_OK|MB_ICONEXCLAMATION "This uninstaller is already running."
      Abort
FunctionEnd

Function un.onUnInstSuccess
   UAC::Unload ;Must call unload!
FunctionEnd

Function un.onUnInstFailed
   UAC::Unload ;Must call unload!
FunctionEnd

# Section "Uninstall"
#   Delete "$INSTDIR\Uninstall.exe"
#   Delete "$INSTDIR\percolator.exe"
#   RMDir "$INSTDIR"
# SectionEnd
