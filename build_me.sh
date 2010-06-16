#!/bin/sh



rm -rf /tmp/build && mkdir /tmp/build && cd /tmp/build  && ~/cmake-2.8.1-Linux-i386/bin/cmake -DCMAKE_INSTALL_PREFIX=/tmp/install '-DCMAKE_PREFIX_PATH=/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu;/scratch/esjolund/xerces-c-3.1.1/inst/;/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu/;/scratch/esjolund/tokyocabinet-1.4.44/inst/' /scratch/e/nypercol/percolator/

export PKG_CONFIG_PATH=/scratch/esjolund/xerces-c-3.1.1/inst/lib/pkgconfig/
export PATH=/scratch/esjolund/xerces-c-3.1.1/inst/bin/:$PATH
export LD_LIBRARY_PATH=/scratch/esjolund/xerces-c-3.1.1/inst/lib/:$LD_LIBRARY_PATH

# as soon as a command fails terminate this script
set -e

percolatorSourceDir=`/bin/pwd`
buildDir=/tmp/percolatorBuildDir
installDir=/tmp/percolatorInstallDir
cd $buildDir 


# you might want to add something like
#'-DCMAKE_PREFIX_PATH=/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu;/scratch/esjolund/xerces-c-3.1.1/inst/;/scratch/esjolund/xsd-3.3.0-x86_64-linux-gnu/;/scratch/esjolund/tokyocabinet-1.4.44/inst/'
# 

cmake -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=FALSE DCMAKE_INSTALL_PREFIX=$installDir  $percolatorSourceDir
make -j8


