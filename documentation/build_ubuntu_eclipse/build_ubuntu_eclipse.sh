#!/bin/sh
# Mattia Tomasoni - Percolator Project
# Script that builds Percolator as an Eclipse Project
# Parameters: none

###############################
# USER MUST SET THESE VARIABLES
###############################
home=/home/mattia # location where the percolator source from the repositories will be stored
buildDir=/home/mattia/percolatorBuild # location where it will be built
installDir=/home/mattia/percolatorInstall # location where it will be installed
codesynthesisDir=/home/mattia/codesynthesis # location for the cosesynthesis library

# install necessary libraries
set -e # as soon as a command fails terminate this script
sudo apt-get install libxerces-c-dev libboost-dev gengetopt libsqlite3-dev cmake libtokyocabinet-dev


# clean old runs of this script
echo ""
echo "**********"
echo "STEP 1: clean old runs..." 
echo "**********"
if [ -d "$codesynthesisDir" ]; then
    rm -r $codesynthesisDir
fi
if [ -d "$buildDir" ]; then
    rm -r $buildDir
fi
if [ -d "$installDir" ]; then
    rm -r $buildDir
fi
if [ -d "$home/percolator" ]; then
    rm -r -f $home/percolator
fi
echo "STEP 1 DONE."
echo ""


# get codesynthesis
echo "**********"
echo "STEP 2: codesynthesis..." 
echo "**********"
mkdir $codesynthesisDir
cd $codesynthesisDir
wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
tar xfj $codesynthesisDir/xsd-3.3.0-x86_64-linux-gnu.tar.bz2
echo "STEP 2 DONE."
echo ""


echo "**********"
echo "STEP 3: get percolator from repositories..." # get percolator
echo "**********"
cd $home
git clone git@github.com:percolator/percolator.git
mkdir $buildDir
echo "STEP 3 DONE."
echo ""


echo "**********"
echo "STEP 4 running cmake..." # run cmake with Eclipse CDT4 option
echo "**********"
cd $buildDir
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug -DGOOGLE_TEST=TRUE -DCMAKE_INSTALL_PREFIX=$installDir -DCMAKE_PREFIX_PATH=$codesynthesisDir/xsd-3.3.0-x86_64-linux-gnu/ $home/percolator
echo "STEP 4 DONE."
echo ""


echo "**********"
echo "SUCCESS!" $buildDir "contains a valid Eclipse Percolator project."
echo "**********"
echo ""
