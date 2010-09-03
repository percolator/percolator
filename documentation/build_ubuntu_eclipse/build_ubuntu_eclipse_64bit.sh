#!/bin/sh
# Mattia Tomasoni - Percolator Project
# Script that builds Percolator as an Eclipse Project
# Parameters: none

src=$(pwd)
buildDir=$src/percolatorBuild # location where it will be built
installDir=$src


# install necessary libraries
set -e # as soon as a command fails terminate this script
sudo apt-get install libxerces-c-dev libboost-dev gengetopt libsqlite3-dev cmake libtokyocabinet-dev git-core zlib1g zlib1g-dev build-essential libgtest-dev


# get codesynthesis
echo "*******"
echo "STEP 1: installing codesynthesis library..."
echo "*******"
cd $src
if [ ! -d "xsd-3.3.0-x86_64-linux-gnu" ]; then
wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
    tar xfj xsd-3.3.0-x86_64-linux-gnu.tar.bz2
    rm xsd-3.3.0-x86_64-linux-gnu.tar.bz2
fi
echo "STEP 1 DONE."
echo ""


echo "*******"
echo "STEP 2: clone/pull percolator from repositories..." # get percolator
echo "*******"
cd $src
if [ ! -d "percolator" ]; then
git clone git@github.com:percolator/percolator.git
else
cd percolator
    git pull
    cd ..
fi
echo "STEP 2 DONE."
echo ""


echo "*******"
echo "STEP 3 running cmake..." # run cmake with Eclipse CDT4 option
echo "*******"
cd $src
echo $src
if [ -d "$buildDir" ]; then
rm -r $buildDir
fi
mkdir $buildDir
cd $buildDir
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug -DGOOGLE_TEST=TRUE -DCMAKE_INSTALL_PREFIX=$installDir -DCMAKE_PREFIX_PATH=$src/xsd-3.3.0-x86_64-linux-gnu/ $src/percolator
echo "STEP 3 DONE."
echo ""


echo "*******"
echo "STEP 4 building, installing and testing..." # run cmake with Eclipse CDT4 option
echo "*******"
cd $buildDir
make
make install
make test
echo "STEP 4 DONE."
echo ""

echo "*******"
echo "SUCCESS!" $buildDir "contains a valid Eclipse project."
echo "*******"
echo ""
