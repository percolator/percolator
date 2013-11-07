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

if [ -z ${build_dir} ]; then
  build_dir="$(mktemp -d --tmpdir ubuntu_build_XXXX)";
fi
if [ -z ${src_dir} ]; then
  if [ -n  ${branch} ]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir ubuntu_build_XXXX)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [ -z ${release_dir} ]; then
  release_dir=${HOME}/release
fi

echo "Building the Percolator packages with src=${src_dir} and build=${build_dir} for the user"
whoami;


# chkconfig sshd on
# usermod lukask -a -G wheel

sudo yum install -y gcc gcc-c++ cmake wget rpm-build
sudo yum install -y tokyocabinet-devel boost boost-devel sqlite-devel zlib-devel 

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

mkdir ${build_dir}
cd ${build_dir}
tar xzf ${src_dir}/${xer}.tar.gz 
cd ${xer}/
#./configure --disable-network --disable-threads --enable-transcoder-gnuiconv --enable-static
./configure --disable-network --disable-threads --enable-static
cd src/
make -j4
ln -s .libs/libxerces-c.a .
ranlib libxerces-c.a

# download, compile and link percolator

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${build_dir}/${xer}/src;${src_dir}/${xsd}/"  ${src_dir}/percolator
make -j4 package
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${build_dir}/${xer}/src;${src_dir}/${xsd}/" ${src_dir}/percolator/src/converters
make -j4 package
make -j4 package

echo "build directory is : "${build_dir}";
cp -v per*.rpm ${release_dir}
