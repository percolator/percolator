#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/mnt/VirtualBoxShare
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild

# the Ubuntu package xsdcxx is as of 2010-10-21 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
# need both version for windows and for linux
cd $base
if [ ! -d "xsd-3.3.0-x86_64-linux-gnu" ]; then
  wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
  tar xfj xsd-3.3.0-x86_64-linux-gnu.tar.bz2
  rm xsd-3.3.0-x86_64-linux-gnu.tar.bz2
fi
if [ ! -d "xsd-3.3.0-i686-windows" ]; then
  wget http://codesynthesis.com/download/xsd/3.3/windows/i686/xsd-3.3.0-i686-windows.zip
  unzip ./xsd-3.3.0-i686-windows.zip
  rm xsd-3.3.0-i686-windows.zip
fi
# the Ubuntu package xsdcxx is as of 2010-10-21 version 1.3 and we want version 1.5
# so we download it from the Googletest home page
if [ ! -d "gtest-1.5.0" ]; then
  wget http://googletest.googlecode.com/files/gtest-1.5.0.tar.gz
  tar xzf gtest-1.5.0.tar.gz
  rm gtest-1.5.0.tar.gz
  cd gtest-1.5.0; mkdir build
  cd build; cmake ../ 
  make
  cd ..
fi

if [ -d "$percolatorBuild" ]; then
  rm -r $percolatorBuild
fi
mkdir $percolatorBuild
cd $percolatorBuild


cmake  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/cmake/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-i686-windows;'$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DEXCLUDE_CONVERTERS=ON -DEXCLUDE_ELUDE=OFF $percolatorSource

make -j 8 win32installer
