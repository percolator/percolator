#!/bin/sh
# as soon as a command fails terminate this script
set -e

percolatorSourceDir=`/bin/pwd`
buildDir=/tmp/percolatorBuildDir
installDir=/tmp/percolatorInstallDir
mkdir -p $buildDir 
cd $buildDir 


# you might want to add something like
#'-DCMAKE_PREFIX_PATH=/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu;/scratch/esjolund/xerces-c-3.1.1/inst/;/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu/;/scratch/esjolund/tokyocabinet-1.4.44/inst/'
# 

cmake -DCMAKE_PREFIX_PATH=~/src/xsd-3.3.0-x86_64-linux-gnu -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=FALSE DCMAKE_INSTALL_PREFIX=$installDir  $percolatorSourceDir

#cmake -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=FALSE DCMAKE_INSTALL_PREFIX=$installDir  $percolatorSourceDir
make -j8


