#!/bin/bash
# Requirements are:
#	Command line tools
#	MacPorts
#	CMake
#----------------------------------------

current_path=$PWD
Ppath="$(dirname ${BASH_SOURCE})"
cd ${Ppath}
Spath=$PWD

if [[ $1 ]]
        then
        release=$1
else 
        release=$HOME/release
fi 

tmp_dir="$(mktemp -d /tmp/OSX_XXXX)"
mkdir ${tmp_dir}/src ${tmp_dir}/build ${tmp_dir}/src/percolator
src="${tmp_dir}/src";build="${tmp_dir}/build";

cp -R "${Spath}"/../../ ${src}/percolator/

sudo port install tokyocabinet boost zlib

#----------------------------------------
xsd=xsd-3.3.0-i686-macosx
cd ${src}
curl -O http://www.codesynthesis.com/download/xsd/3.3/macosx/i686/${xsd}.tar.bz2
tar -xjf ${xsd}.tar.bz2
sed -i -e 's/setg/this->setg/g' ${xsd}/libxsd/xsd/cxx/zc-istream.txx
sed -i -e 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/parsing.txx
sed -i -e 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx

extr=${xsd}/libxsd/xsd/cxx/tree/xdr-stream-extraction.hxx
inse=${xsd}/libxsd/xsd/cxx/tree/xdr-stream-insertion.hxx
sed -i -e 's/ uint8_t/ unsigned char/g' ${extr} ${inse}
sed -i -e 's/ int8_t/ char/g' ${extr} ${inse}
sed -i -e 's/xdr_int8_t/xdr_char/g' ${extr} ${inse}
sed -i -e 's/xdr_uint8_t/xdr_u_char/g' ${extr} ${inse}
#------------------------------------------
xer=xerces-c-3.1.1
mkdir -p ${build}
cd ${build}
if [[ -d /usr/local/include/xercesc ]]
	then
	echo "Xerces is already installed."
else
        cd ${src}
	curl -O http://apache.mirrors.spacedump.net//xerces/c/3/sources/${xer}.tar.gz
        cd ${build}
	tar xzf ${src}/${xer}.tar.gz 
	cd ${xer}/
	./configure CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64" --disable-network --disable-threads
	make
	sudo make install
fi

#-------------------------------------------
mkdir -p ${build}/percolator
cd ${build}/percolator

cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${src}/${xsd}/"  ${src}/percolator
make package

mkdir -p ${build}/converters
cd ${build}/converters

cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${src}/${xsd}/"  -DSERIALIZE="TokyoCabinet" ${src}/percolator/src/converters
make package
#--------------------------------------------
mkdir -p ${release}
cp ${build}/percolator/*.dmg ${release}
cp ${build}/converters/*.dmg ${release}
