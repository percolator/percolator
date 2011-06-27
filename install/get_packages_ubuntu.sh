#!/bin/sh
# as soon as a command fails terminate this script

base=/tmp

#make sure user has sudo rights
if [ "$(id -u)" != "0" ]; then
   echo "This script must be run as root" 1>&2
   exit 1
fi
set -e

#install libraries
apt-get install libxerces-c-dev libboost-dev build-essential cmake

#detect kernel architecture
if [ "$(uname -a | grep x86_64)" != "" ]; then
  architecture="64"
else
  architecture="32"
fi

#download xsd and googletest
cs_version="3.3"
cs_site="http://codesynthesis.com/download/xsd/"${cs_version}"/linux-gnu"
if [ ${architecture} = "64" ]; then
  cs_site=${cs_site}/"x86_64"
  cs_pack="xsd-"${cs_version}".0-x86_64-linux-gnu"
else
  cs_site=${cs_site}/"i686"
  cs_pack="xsd-"${cs_version}".0-i686-linux-gnu"
fi

gt_site="http://googletest.googlecode.com/files"
gt_pack="gtest-1.5.0"

my_wd=`pwd`
cd $base

#This is a temporal thing Codesynthesis 3.3 is to be part of the next ubuntu release
echo "Fetching codesynthesis, storing it in ${base}"
wget -nc "${cs_site}/${cs_pack}.tar.bz2"
tar xjf "${cs_pack}.tar.bz2"
rm "${cs_pack}.tar.bz2"

echo "Fetching GoogleTest, storing it in ${base}"
wget -nc "${gt_site}/${gt_pack}.tar.gz"
tar xzf "${gt_pack}.tar.gz"
rm "${gt_pack}.tar.gz"

echo "Building GoogleTest"
cd gtest-1.5.0; mkdir -p build
cd build; cmake ../ 
make
cd ..

echo "To build, invoke cmake from the build directory with the following option:\n"
echo "-DCMAKE_PREFIX_PATH=\"${base}/${cs_pack};${base}/${gt_pack}\""
