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
  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir build_XXXX)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} and build=${build_dir} for the user"
whoami;

sudo dnf -y install gcc gcc-c++ wget rpm-build cmake
sudo dnf install sqlite-devel zlib-devel bzip2-devel
sudo dnf install tokyocabinet-devel xerces-c-devel
sudo dnf install libtirpc libtirpc-devel
sudo dnf -y --enablerepo=powertools install gtest

cd ${src_dir}

# read all urls and file names from a centralized kb file
source percolator/admin/builders/_urls_and_file_names_.sh

cd ${build_dir}

# download and install boost for CentOS < 8, since these install boost <= 1.56 which has problems with the header only library includes
if [[ $(rpm -q --queryformat '%{VERSION}' centos-release) < 8 ]]; then
  if [ ! -d ${centos_boost} ]; then
    echo "  Installing boost"
    wget --quiet -O ${centos_boost}.tar.bz2 ${centos_boost_url}
    tar xjf ${centos_boost}.tar.bz2
    cd ${centos_boost}/
    ./bootstrap.sh
    ./b2 address-model=64 threading=multi -j4 --with-system --with-filesystem --with-serialization -d0
  fi
else
  sudo yum install -y boost-static boost-devel
fi

# download and install xsd

if hash xsdcxx 2>/dev/null; then
  echo "  XSD has been installed previously, remove if you want a clean install"
else
  echo "  Installing XSD"
  wget --quiet ${centos_xsd_url}
  rpm -ivh ${centos_xsd}.rpm
fi

# download, compile and link percolator

echo "Installing percolator"

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${build_dir}/${centos_boost}" -DXML_SUPPORT=OFF ${src_dir}/percolator
make -j 4;
make -j 4 package;

# Fix to handle alt. rpc location
# export CFLAGS=`pkg-config --cflags libtirpc`
# export CXXFLAGS=-I/usr/include/tirpc

mkdir -p ${build_dir}/percolator;
cd ${build_dir}/percolator;
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${build_dir}/${centos_boost}" -DXML_SUPPORT=ON ${src_dir}/percolator;
make -j 4;
make -j 4 package;

#-----cmake-----
mkdir -p ${build_dir}/percolator-debug;
cd $build_dir/percolator-debug;
cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Debug -DGOOGLE_TEST=1 -DCMAKE_INSTALL_PREFIX=./local-usr -DCMAKE_PREFIX_PATH="${build_dir}/${centos_boost}" -DXML_SUPPORT=ON ${src_dir}/percolator;
#-----make------
make -j 4;


mkdir -p ${build_dir}/converters
cd ${build_dir}/converters
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${build_dir}/${centos_boost}" -DSERIALIZE="TokyoCabinet" ${src_dir}/percolator/src/converters
make -j 4;
make -j 4 package;

mkdir -p ${build_dir}/elude
cd ${build_dir}/elude
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="${build_dir}/${centos_boost}" ${src_dir}/percolator/src/elude_tool
make -j 4;
make -j 4 package;

echo "build directory was : ${build_dir}";

cp -v ${build_dir}/{percolator-noxml,percolator,converters,elude}/*.rpm ${release_dir};
