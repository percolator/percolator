#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/tmp
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild
percolatorInstall=$base/percolatorInstall

# the Ubuntu package xsdcxx is as of 2010-10-21 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
cd $base
if [ ! -d "xsd-3.3.0-x86_64-linux-gnu" ]; then
  wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
  tar xfj xsd-3.3.0-x86_64-linux-gnu.tar.bz2
  rm xsd-3.3.0-x86_64-linux-gnu.tar.bz2
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
if [ -d "$percolatorInstall" ]; then
rm -r $percolatorInstall
fi
mkdir $percolatorInstall
cd $percolatorBuild

cmake -G"Eclipse CDT4 - Unix Makefiles" -DSTATIC=off -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=TRUE -DEXCLUDE_CONVERTERS=FALSE -DEXCLUDE_ELUDE=FALSE -DCMAKE_INSTALL_PREFIX=$percolatorInstall -DCMAKE_PREFIX_PATH=/tmp/xsd-3.3.0-x86_64-linux-gnu/ $percolatorSource

make -j 8
make install

PATH=$PATH:$percolatorInstall
export PATH

make test
