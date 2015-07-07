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


# chkconfig sshd on
# usermod lukask -a -G wheel

sudo yum install -y gcc gcc-c++ wget rpm-build
sudo yum install -y sqlite-devel zlib-devel bzip2-devel
# sudo yum install -y tokyocabinet-devel boost-static # are not available for CentOS5.10
# sudo yum install -y boost-devel cmake # outdated versions in CentOS5.10

cd ${src_dir}

# download and install cmake
cm=cmake-2.8.12.2-Linux-i386
cmake_bin=""
if [ -e ${cm}.sh ]; then
  echo "  CMake has been installed previously"
else
  echo "  Installing CMake"
  wget --quiet http://www.cmake.org/files/v2.8/${cm}.sh
  mkdir ${cm}
  sudo sh ./${cm}.sh --skip-license --prefix=${cm}
fi
cmake_bin="${src_dir}/${cm}/bin/cmake"

# download, compile and link xerces
xer=xerces-c-3.1.1
if [ -e ${xer}.tar.gz ]; then
  echo "  XercesC has been installed previously, remove if you want a clean install"
else
  echo "  Installing XercesC"
  wget --quiet http://archive.apache.org/dist/xerces/c/3/sources/${xer}.tar.gz
  tar xzf ${xer}.tar.gz 
  cd ${xer}
  ./configure --disable-network --disable-threads --disable-shared --enable-static
  cd src
  make -j 4
  ln -s .libs/libxerces-c.a .
  ranlib libxerces-c.a
  cd ../../
fi

# download and patch xsd
xsd=xsd-4.0.0+dep
if [ -e ${xsd}.tar.bz2 ]; then
  echo "  XSD has been installed previously, remove if you want a clean install"
else
  echo "  Installing XSD"
  wget --quiet http://www.codesynthesis.com/download/xsd/4.0/${xsd}.tar.bz2
  tar xjf ${xsd}.tar.bz2
  cd ${xsd}
  make CPPFLAGS=-I../${xer}/src LDFLAGS=-L../${xer}/src/.libs
  ./xsd/xsd/xsd --version
  mkdir -p xsd/bin
  mv ./xsd/xsd/xsd xsd/bin
  cd ../
  #sed -i 's/setg/this->setg/g' ${xsd}/libxsd/xsd/cxx/zc-istream.txx
  #sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/parsing.txx
  #sed -i 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx
fi

# download, compile and link tokyocabinet
tokyo=tokyocabinet-1.4.48
if [ -e ${tokyo}.tar.gz ]; then
  echo "  TokyoCabinet has been installed previously, remove if you want a clean install"
else
  echo "  Installing TokyoCabinet"
  wget --quiet http://fallabs.com/tokyocabinet/${tokyo}.tar.gz
  tar zxf ${tokyo}.tar.gz
  cd ${tokyo}
  ./configure
  make -j 4
  make install
  cd ../
fi

boost=boost_1_57_0
boost_version="1.57.0"
if [ -e ${boost}.tar.gz ]; then
  echo "  Boost ${boost_version} has been installed previously, remove if you want a clean install"
else
  echo "  Installing Boost, this might take a while..."
  wget http://sourceforge.net/projects/boost/files/boost/${boost_version}/${boost}.tar.gz/download
  tar zxf ${boost}.tar.gz
  cd ${boost}
  sudo sh ./bootstrap.sh
  sudo ./bjam --layout=system threading=multi --with-system --with-filesystem --with-serialization -d1 install
  cd ../
fi

# download, compile and link percolator

echo "Installing percolator"

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml
${cmake_bin} -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF ${src_dir}/percolator
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator
${cmake_bin} -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=ON -DBOOST_ROOT=/usr/local -DBOOST_LIBRARYDIR=/usr/local/lib -DCMAKE_PREFIX_PATH="${src_dir}/${xer}/src;${src_dir}/${xsd}/xsd/" ${src_dir}/percolator
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters
${cmake_bin} -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DBOOST_ROOT=/usr/local -DBOOST_LIBRARYDIR=/usr/local/lib -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${src_dir}/${xer}/src;${src_dir}/${xsd}/xsd/" ${src_dir}/percolator/src/converters
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/elude
cd ${build_dir}/elude
${cmake_bin} -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr ${src_dir}/percolator/src/elude_tool
make -j 4;
make -j 4 package;
cp elude*.rpm ${release_dir}

echo "build directory was : ${build_dir}";

