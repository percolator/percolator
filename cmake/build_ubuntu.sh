#!/bin/sh
# as soon as a command fails terminate this script
set -e

base=/tmp
percolatorSource=`/bin/pwd`
percolatorBuild=$base/percolatorBuild

if [ -d "$percolatorBuild" ]; then
  rm -r $percolatorBuild
fi
mkdir $percolatorBuild
cd $percolatorBuild

cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX=$base $percolatorSource
make -j8
make install
PATH=$PATH:$base/bin; export PATH
make test
