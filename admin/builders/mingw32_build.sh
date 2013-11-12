#!/bin/bash

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
  build_dir="$(HOME)/build";
#  build_dir="$(mktemp -d --tmpdir build_XXXX)";
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


sudo yum install -y cmake wget mingw-w64-tools mingw32-gcc-c++ mingw32-filesystem mingw-binutils-generic mingw32-nsis
sudo yum install -y mingw32-boost-static mingw32-sqlite mingw32-zlib mingw32-curl mingw32-pthreads

cd ${src_dir}

# download and patch xsd

xsd=xsd-3.3.0-x86_64-linux-gnu
wget --quiet http://www.codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/${xsd}.tar.bz2
tar xjf ${xsd}.tar.bz2
sed -i 's/setg/this->setg/g' ${xsd}/libxsd/xsd/cxx/zc-istream.txx
sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/parsing.txx
sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx

# download, compile and link xerces
xer=xerces-c-3.1.1

wget --quiet http://apache.mirrors.spacedump.net//xerces/c/3/sources/${xer}.tar.gz

mkdir -p ${build_dir}
cd ${build_dir}

tar xzf ${src_dir}/${xer}.tar.gz 
cd ${xer}/
./configure --disable-network --disable-threads --enable-transcoder-windows --disable-static --enable-shared --host=i686-w64-mingw32 --prefix=/usr/i686-w64-mingw32/sys-root/mingw
cd src/
make libxerces_c_la_LDFLAGS="-release 3.1 -no-undefined" -j4
sudo make install

# download, compile and link percolator

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator


mingw32-cmake -DCMAKE_PREFIX_PATH="${src_dir}/${xsd}/;${src_dir}/${xer}/src/" -DCMAKE_BUILD_TYPE=Release ${src_dir}/percolator
make -j4 package

cp -v per*.exe ${release_dir}
 
mkdir -p ${build_dir}/converters
cd ${build_dir}/converters

mingw32-cmake -DSERIALIZE="Boost" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${src_dir}/${xsd}/" ${src_dir}/percolator/src/converters
make -j4 package
cp -v per*.exe ${release_dir}

echo "build directory is : ${build_dir}";
