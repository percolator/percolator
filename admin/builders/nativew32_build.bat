call "C:\Program Files\Microsoft Visual Studio 12.0\Common7\Tools\VsDevCmd.bat"
call net use z: "\\VBOXSVR\VagrantWin7"
call cd /D Z:\percolator
call md build
call cd build
call cmake -G "Visual Studio 12" ..
call msbuild PERCOLATOR.sln /p:VCTargetsPath="C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V120"
