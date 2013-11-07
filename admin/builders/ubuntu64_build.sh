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

#------------------------------------------------------------------------
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

sudo apt-get update;
sudo apt-get upgrade;
#sudo apt-get -y install g++ make cmake rpm fakeroot;
sudo apt-get -y install g++ make rpm fakeroot;
# Need a never copy of cmake
# Remove the secion below once they updated cmake
cd ${src_dir}
wget -q http://www.cmake.org/files/v2.8/cmake-2.8.12.tar.gz
tar xzf cmake-2.8.12.tar.gz
cd cmake-2.8.12/
./bootstrap;
make -j2; 
sudo make install;
# end of section to remove
sudo apt-get -y install xsdcxx libxerces-c-dev libboost-dev libboost-filesystem-dev;
sudo apt-get -y install libboost-system-dev libboost-thread-dev libsqlite3-dev libtokyocabinet-dev zlib1g-dev;

#------------------------------------------------------------------------
mkdir -p $build_dir/percolator $build_dir/converters;

######percolator########
#-----cmake-----
cd $build_dir/percolator;
echo -n "cmake percolator.....";
fakeroot -- cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr $src_dir/percolator;
#-----make------
echo -n "make percolator (this will take few minutes).....";
fakeroot -- make -j2;
make -j2 package;

#######converters########
cd $build_dir/converters
#-----cmake-----
echo -n "cmake converters.....";
fakeroot -- cmake -DTARGET_ARCH=amd64 -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DSERIALIZE="TokyoCabinet" $src_dir/percolator/src/converters;

#-----make------
echo -n "make converters (this will take few minutes).....";

fakeroot -- make -j2;
make -j2 package;

###########################
cp $build_dir/{percolator,converters}/*.deb ${release_dir};
echo "Finished buildscript execution";
echo "in build directory ${build_dir}";
