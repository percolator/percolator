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
cs_version="3.3"
cs_site="http://codesynthesis.com/download/xsd/"${cs_version}"/linux-gnu"
if [ ! -d "xsd-3.3.0-x86_64-linux-gnu" ]; then
  cs_site64=${cs_site}/"x86_64"
  cs_pack="xsd-"${cs_version}".0-x86_64-linux-gnu"
  wget -nc "${cs_site64}/${cs_pack}.tar.bz2"
  tar xjf "${cs_pack}.tar.bz2"
  rm "${cs_pack}.tar.bz2"
fi
if [ ! -d "xsd-3.3.0-i686-windows" ]; then
  cs_site32=${cs_site}/"i686"
  cs_pack="xsd-"${cs_version}".0-i686-linux-gnu"
  wget -nc "${cs_site32}/${cs_pack}.tar.bz2"
  tar xjf "${cs_pack}.tar.bz2"
  rm "${cs_pack}.tar.bz2"
fi

if [ -d "$percolatorBuild" ]; then
  rm -r $percolatorBuild
fi
mkdir $percolatorBuild
cd $percolatorBuild


cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/cmake/windows_percolator/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-i686-windows;'$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $percolatorSource

make -j 8 win32installer
