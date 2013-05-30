#!/bin/bash

post="_mingw64"
branch="branch-2-05"

sudo yum install -y cmake wget mingw-w64-tools mingw-binutils-generic
sudo yum install -y mingw64-tokyocabinet mingw64-boost mingw64-sqlite mingw64-zlib

src=/tmp/src${post}
build=/tmp/build${post}

rm -fr ${src} ${build}
mkdir -p ${src} ${build}


cd ${src}

# download and patch xsd

xsd=xsd-3.3.0-x86_64-linux-gnu
wget http://www.codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/${xsd}.tar.bz2
tar xvjf ${xsd}.tar.bz2
sed -i 's/setg/this->setg/g' ${xsd}/libxsd/xsd/cxx/zc-istream.txx
sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/parsing.txx
sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx

# download, compile and link xerces
xer=xerces-c-3.1.1

wget http://apache.mirrors.spacedump.net//xerces/c/3/sources/${xer}.tar.gz

cd ${build}

tar xvzf ${src}/${xer}.tar.gz 
cd ${xer}/
./configure --disable-network --disable-threads --enable-transcoder-windows --en
able-shared --host=x86_64-w64-mingw32
cd src/
make libxerces_c_la_LDFLAGS="-release 3.1 -no-undefined" -j4
sudo make install

# download, compile and link percolator

cd  ${src}
git clone git://github.com/percolator/percolator.git
cd percolator
git checkout ${branch}

mkdir -p ${build}/percolator
cd ${build}/percolator


mingw64-cmake -DCMAKE_PREFIX_PATH="${src}/${xsd}/"  ${src}/percolator
make -j4 package
 

mkdir -p ${build}/converters
cd ${build}/converters

mingw64-cmake -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="${src}/${xsd}/" ${src}/percolator/src/converters
make -j4 package
