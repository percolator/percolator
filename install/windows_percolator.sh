#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/tmp
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild
percolatorInstall=$base/percolatorInstall

# the Ubuntu package xsdcxx is as of 2010-06-16 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
# need both version for windows and for linux
cd $base
wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
tar xfj ./xsd-3.3.0-x86_64-linux-gnu.tar.bz2
rm ./xsd-3.3.0-x86_64-linux-gnu.tar.bz2
wget http://codesynthesis.com/download/xsd/3.3/windows/i686/xsd-3.3.0-i686-windows.zip
unzip ./xsd-3.3.0-i686-windows.zip
rm ./xsd-3.3.0-i686-windows.zip
mkdir $percolatorBuild
mkdir $percolatorInstall
cd $percolatorBuild

cmake -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/install/cmake/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-i686-windows;'$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DEXCLUDE_CONVERTERS=ON $percolatorSource

make -j 8 win32installer
make install

PATH=$PATH:$percolatorInstall
export PATH

make test
