#!/bin/bash

# managing input arguments

if [ $# -eq 2 ];
  then src=$1;build=$2;mkdir -p ${build};
elif [ $# -eq 1 ];
  then sudo yum git;
  tmp_dir="$(mktemp -d --tmpdir mingw_tmp_XXXX)";
  mkdir ${tmp_dir}/src;mkdir ${tmp_dir}/build;
  src="${tmp_dir}/src";build="${tmp_dir}/build";
  git clone --branch "$1" https://github.com/percolator/percolator.git ${src}/percolator;
else 
  echo "Please add either one argument as branch name or two arguments for your source directory containing percolator/ and build directory";
  exit 1;
fi;


sudo yum install -y cmake wget mingw-w64-tools mingw64-filesystem mingw-binutils-generic mingw32-nsis
sudo yum install -y mingw64-boost mingw64-sqlite mingw64-zlib mingw64-curl mingw64-pthreads

#release=${home}/rel

cd ${src}

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

mkdir -p ${build}
cd ${build}

tar xzf ${src}/${xer}.tar.gz 
cd ${xer}/
./configure --disable-network --disable-threads --enable-shared --host=x86_64-w64-mingw32 --prefix=/usr/x86_64-w64-mingw32/sys-root/mingw
#./configure --disable-network --disable-threads --enable-transcoder-windows --enable-shared --host=x86_64-w64-mingw32 --prefix=/usr/x86_64-w64-mingw32/sys-root/mingw
cd src/
make libxerces_c_la_LDFLAGS="-release 3.1 -no-undefined" -j4
sudo make install

# download, compile and link percolator

mkdir -p ${build}/percolator
cd ${build}/percolator


mingw64-cmake -DCMAKE_PREFIX_PATH="${src}/${xsd}/;${src}/${xer}/src/"  ${src}/percolator
make -j4 package

cp per*.exe ${rel}
 

mkdir -p ${build}/converters
cd ${build}/converters

mingw64-cmake -DSERIALIZE="Boost" -DCMAKE_PREFIX_PATH="${src}/${xsd}/" ${src}/percolator/src/converters
make -j4 package

echo "build directory is : ${build}";
cp per*.exe ${rel}
