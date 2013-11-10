#!/bin/bash
# Requirements are:
#	Command line tools
#	MacPorts
#	CMake
#----------------------------------------


# managing input arguments
while getopts “s:b:r:t:” OPTION; do
  case $OPTION in
    s) src_dir=${OPTARG};;
    t) branch=${OPTARG};;
    r) release_dir=${OPTARG};;
    b) build_dir=${OPTARG};;
    \?) echo "Invalid option: -${OPTARG}" >&2;;
  esac
done

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir src_XXXX)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for the user"
whoami

sudo port install tokyocabinet boost zlib

#----------------------------------------
xsd=xsd-3.3.0-i686-macosx
cd ${src_dir}
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
mkdir -p ${build_dir}
cd ${build_dir}
if [[ -d /usr/local/include/xercesc ]]
	then
	echo "Xerces is already installed."
else
    cd ${src_dir}
	curl -O http://apache.mirrors.spacedump.net/xerces/c/3/sources/${xer}.tar.gz
    cd ${build_dir}
	tar xzf ${src_dir}/${xer}.tar.gz
	cd ${xer}/
	./configure CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64" --disable-network --disable-threads
	make
	sudo make install
fi

#-------------------------------------------
mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${src_dir}/${xsd}/"  ${src_dir}/percolator
make -j 2
make -j 2 package

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters

cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${src_dir}/${xsd}/"  -DSERIALIZE="TokyoCabinet" ${src_dir}/percolator/src/converters
make -j 2
make -j 2 package
#--------------------------------------------
mkdir -p ${release_dir}
cp -v ${build_dir}/percolator/*.dmg ${release_dir}
cp -v ${build_dir}/converters/*.dmg ${release_dir}
