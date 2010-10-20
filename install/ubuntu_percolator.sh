#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/tmp
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild
percolatorInstall=$base/percolatorInstall

# the Ubuntu package xsdcxx is as of 2010-06-16 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
cd $base
wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
tar xfj ./xsd-3.3.0-x86_64-linux-gnu.tar.bz2
rm ./xsd-3.3.0-x86_64-linux-gnu.tar.bz2
mkdir $percolatorBuild
mkdir $percolatorInstall
cd $percolatorBuild

cmake -G"Eclipse CDT4 - Unix Makefiles" -DSTATIC=off -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=FALSE -DEXCLUDE_CONVERTERS=TRUE -DEXCLUDE_ELUDE=TRUE -DCMAKE_INSTALL_PREFIX=$percolatorInstall -DCMAKE_PREFIX_PATH=/tmp/xsd-3.3.0-x86_64-linux-gnu/ $percolatorSource

make -j 8
make install

PATH=$PATH:$percolatorInstall
export PATH

make test
