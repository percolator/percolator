#!/bin/sh
# as soon as a command fails terminate this script

base=/scratch

#detect (kernel) architecture
if [ "$(uname -a | grep x86_64)" != "" ]; then
  architecture="64"
else
  architecture="32"
fi

if [ ${architecture} = "64" ]; then
  cs_site="http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64"
  cs_pack="xsd-3.3.0-x86_64-linux-gnu"
else
  cs_site="http://codesynthesis.com/download/xsd/3.3/linux-gnu/i686"
  cs_pack="xsd-3.3.0-i686-linux-gnu"
fi

gt_site="http://googletest.googlecode.com/files"
gt_pack="gtest-1.5.0"

if [ "$(id -u)" != "0" ]; then
   echo "This script must be run as root" 1>&2
   exit 1
fi

set -e

apt-get install libxerces-c-dev libboost-dev build-essential cmake
apt-get install gengetopt libtokyocabinet-dev zlib1g-dev

my_wd=`pwd`
cd $base

#This is a temporal thing Codesynthesis 3.3 is to be part of the next ubuntu release
echo "Fetching codesynthesis, storing it in ${base}"

wget "${cs_site}/${cs_pack}.tar.bz2"
tar xjf "${cs_pack}.tar.bz2"
rm "${cs_pack}.tar.bz2"

echo "Fetching GoogleTest, storing it in ${base}"

wget "${gt_site}/${gt_pack}.tar.gz"
tar xzf "${gt_pack}.tar.gz"
rm "${gt_pack}.tar.gz"

echo "Building GoogleTest"
cd gtest-1.5.0; mkdir -p build
cd build; cmake ../ 
make
cd ..

echo "To build, invoke cmake from the build directory with the following option:\n"
echo "-DCMAKE_PREFIX_PATH=\"${base}/${cs_pack};${base}/${gt_pack}\""
