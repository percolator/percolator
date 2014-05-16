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

sudo yum install -y gcc gcc-c++ cmake wget rpm-build
sudo yum install -y tokyocabinet-devel boost-static boost-devel sqlite-devel zlib-devel bzip2-devel

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
make -j 4
ln -s .libs/libxerces-c.a .
ranlib libxerces-c.a

# download, compile and link percolator

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF ${src_dir}/percolator
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=ON -DCMAKE_PREFIX_PATH="${build_dir}/${xer}/src;${src_dir}/${xsd}/"  ${src_dir}/percolator
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" -DCMAKE_PREFIX_PATH="${build_dir}/${xer}/src;${src_dir}/${xsd}/" ${src_dir}/percolator/src/converters
make -j 4;
make -j 4 package;
cp per*.rpm ${release_dir}

mkdir -p ${build_dir}/elude
cd ${build_dir}/elude
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr ${src_dir}/percolator/src/elude_tool
make -j 4;
make -j 4 package;
cp elude*.rpm ${release_dir}

echo "build directory was : ${build_dir}";


# Fixes permission denied /usr folders bug in CMake 2.8.11 on Fedora 18
# This bug was fixed from version 2.8.12 on, so no need to run it there.
cd ${release_dir}
if [[ `cmake --version | cut -f3 -d' ' | awk -F'.' '{print $2*100+$3}'` -le "811" ]]; then
  wget -q -O rpmrebuild.rpm http://sourceforge.net/projects/rpmrebuild/files/rpmrebuild/2.11/rpmrebuild-2.11-1.noarch.rpm/download 
  sudo rpm -i rpmrebuild.rpm
  rm rpmrebuild.rpm

  backupdir=original_packages
  mkdir $backupdir
  mv *.rpm $backupdir
  cd $backupdir
  for rpm in *.rpm
  do
    rpmrebuild --package --change-spec-whole="sed '/%dir %attr(0755, root, root) \"\/usr\"/d' | sed '/%dir %attr(0755, root, root) \"\/usr\/bin\"/d' | sed '/%dir %attr(0755, root, root) \"\/usr\/share\"/d'" --directory=${release_dir} $rpm
  done
  
  cd ${release_dir}
  mv ./x86_64/*.rpm ./
  rm -rf x86_64
fi

