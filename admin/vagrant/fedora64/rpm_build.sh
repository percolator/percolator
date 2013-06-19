#!/bin/bash

post="Fedora"
src="/vagrant/src"
build="/vagrant/build"


# chkconfig sshd on
# usermod lukask -a -G wheel

yum install -y gcc 
yum install -y gcc-c++ 
yum install -y cmake wget rpm-build
yum install -y tokyocabinet-devel boost boost-devel sqlite-devel zlib-devel 

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

mkdir ${build}
cd ${build}
tar xvzf ${src}/${xer}.tar.gz 
cd ${xer}/
./configure --disable-network --disable-threads --enable-transcoder-gnuiconv --enable-static
cd src/
make -j4
ln -s .libs/libxerces-c.a .
ranlib libxerces-c.a

# download, compile and link percolator

mkdir -p ${build}/percolator
cd ${build}/percolator

cmake -DTARGET_ARCH=amd64 -DCMAKE_PREFIX_PATH="${build}/${xer}/src;${src}/${xsd}/"  ${src}/percolator
make -j4 package
#cp per*.rpm ${rel}

mkdir -p ${build}/converters
cd ${build}/converters
cmake -DTARGET_ARCH=amd64 -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${build}/${xer}/src;${src}/${xsd}/" ${src}/percolator/src/converters
make -j4 package
make -j4 package

#cp per*.rpm ${rel}
