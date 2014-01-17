call "C:\Program Files\Microsoft Visual Studio 12.0\Common7\Tools\VsDevCmd.bat"
call net use z: "\\VBOXSVR\VagrantWin7"
call cd /D Z:\percolator
call md build
call cd build
call cmake -G "Visual Studio 12" -DBOOST_ROOT="C:\\local\\boost_1_55_0" -DBOOST_LIBRARYDIR="C:\\local\\boost_1_55_0\\lib32-msvc-12.0" -DCMAKE_PREFIX_PATH="C:\\local\\xerces-c-3.1.1-x86-windows-vc-10.0;C:\\local\\xsd-3.3.0-i686-windows" -DXML_SUPPORT=ON ..
::call msbuild PERCOLATOR.sln /p:VCTargetsPath="C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V120" /p:Configuration=Release
call msbuild INSTALL.vcxproj /p:VCTargetsPath="C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V120" /p:Configuration=Release
call cd ..\admin\builders
