#!/bin/bash

# managing input arguments


if [ $# -eq 2 ];
  then src=$1;build=$2;
elif [ $# -eq 1 ];
  then sudo yum install git;
  tmp_dir="$(mktemp -d --tmpdir rpm_tmp_XXXX)";
  mkdir ${tmp_dir}/src;mkdir ${tmp_dir}/build;
  src="${tmp_dir}/src";build="${tmp_dir}/build";
  git clone --branch "$1" https://github.com/percolator/percolator.git ${src}/percolator;
else 
  echo "Please add either one argument as branch name or two arguments for your source directory containing percolator/ and build directory";
  exit 1;
fi;

echo "Building the Percolator packages with src=${src} and build=${build}"

# chkconfig sshd on
# usermod lukask -a -G wheel

sudo yum install -y gcc gcc-c++ cmake wget rpm-build
sudo yum install -y tokyocabinet-devel boost boost-devel sqlite-devel zlib-devel 

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

mkdir ${build}
cd ${build}
tar xzf ${src}/${xer}.tar.gz 
cd ${xer}/
#./configure --disable-network --disable-threads --enable-transcoder-gnuiconv --enable-static
./configure --disable-network --disable-threads --enable-static
cd src/
make -j4
ln -s .libs/libxerces-c.a .
ranlib libxerces-c.a

# download, compile and link percolator

mkdir -p ${build}/percolator
cd ${build}/percolator

cmake -DTARGET_ARCH=x86_64 -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${build}/${xer}/src;${src}/${xsd}/"  ${src}/percolator
make -j4 package
#cp per*.rpm ${rel}

mkdir -p ${build}/converters
cd ${build}/converters
cmake -DTARGET_ARCH=x86_64 -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${build}/${xer}/src;${src}/${xsd}/" ${src}/percolator/src/converters
make -j4 package
make -j4 package

echo "build directory is : "${build}";
#cp per*.rpm ${rel}
