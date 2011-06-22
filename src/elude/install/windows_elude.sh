#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/mnt/VirtualBoxShare
eludeSource=`/bin/pwd`
eludeBuild=$base/eludeBuild

# so we download from the Googletest home page
cd $base
if [ ! -d "gtest-1.5.0" ]; then
  wget http://googletest.googlecode.com/files/gtest-1.5.0.tar.gz
  tar xzf gtest-1.5.0.tar.gz
  rm gtest-1.5.0.tar.gz
  cd gtest-1.5.0; mkdir build
  cd build; cmake ../ 
  make
  cd ..
fi

if [ -d "$eludeBuild" ]; then
  rm -r $eludeBuild
fi
mkdir $eludeBuild
cd $eludeBuild


cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$eludeSource/cmake/windows_elude/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $eludeSource

make -j 8 win32installer
