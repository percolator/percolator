#!/bin/bash

post="_mingw64"
src="/vagrant/src"
build="/vagrant/build"


yum install -y cmake wget mingw-w64-tools mingw64-filesystem mingw-binutils-generic mingw32-nsis
yum install -y mingw64-boost mingw64-sqlite mingw64-zlib mingw64-curl mingw64-pthreads

#release=${home}/rel

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
./configure --disable-network --disable-threads --enable-transcoder-windows --enable-shared --host=x86_64-w64-mingw32 --prefix=/usr/x86_64-w64-mingw32/sys-root/mingw
cd src/
make libxerces_c_la_LDFLAGS="-release 3.1 -no-undefined" -j4
make install

# download, compile and link percolator

mkdir -p ${build}/percolator
cd ${build}/percolator


mingw64-cmake -DCMAKE_PREFIX_PATH="${src}/${xsd}/;${src}/${xer}/src/"  ${src}/percolator
make -j4 package

#cp per*.exe ${rel}
 

mkdir -p ${build}/converters
cd ${build}/converters

mingw64-cmake -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="${src}/${xsd}/" ${src}/percolator/src/converters
make -j4 package

#cp per*.exe ${rel}
