#!/bin/sh
# as soon as a command fails terminate this script
set -e

#base=/mnt/VirtualBoxShare
base=/tmp/test
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild

download(){
  site=$1
  pack=$2
  wget -nc "${site}/${pack}.tar.bz2"
  tar xjf "${pack}.tar.bz2"
  rm "${pack}.tar.bz2"
}

# the Ubuntu package xsdcxx is as of 2010-10-21 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
# need both version for windows and for linux
cd $base
cs_version="3.3"
cs_site="http://codesynthesis.com/download/xsd/"${cs_version}"/linux-gnu"
if [ ! -d "xsd-3.3.0-x86_64-linux-gnu" ]; then
  cs_pack="xsd-"${cs_version}".0-x86_64-linux-gnu"
  download $cs_site/"x86_64" $cs_pack
fi
if [ ! -d "xsd-3.3.0-i686-windows" ]; then
  cs_pack="xsd-"${cs_version}".0-i686-linux-gnu"
  download $cs_site/"i686" $cs_pack
fi

if [ -d "$percolatorBuild" ]; then
  rm -r $percolatorBuild
fi
mkdir $percolatorBuild
cd $percolatorBuild


cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/cmake/windows_percolator/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-i686-windows;'$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $percolatorSource

make -j 8 win32installer
