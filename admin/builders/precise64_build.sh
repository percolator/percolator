#!/bin/bash
# managing input arguments

if [ $# -eq 2 ];
  then src_dir=$1;build_dir=$2;
elif [ $# -eq 1 ];
  then sudo apt-get install git;
  tmp_dir="$(mktemp -d --tmpdir precise_tmp_XXXX)";
  mkdir ${tmp_dir}/src;mkdir ${tmp_dir}/build;
  src_dir="${tmp_dir}/src";build_dir="${tmp_dir}/build";
  git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
else 
  echo "Please add either one argument as branch name or two arguments for your source directory containing percolator/ and build directory";
  exit 1;
fi;

echo "Building the Percolator packages with src=${src_dir} and build=${build_dir}"


#------------------------------------------------------------------------
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

sudo apt-get update;
sudo apt-get -y install g++ make cmake rpm fakeroot;
sudo apt-get -y install xsdcxx libxerces-c-dev libboost-dev libboost-filesystem-dev;
#sudo apt-get -y install libboost-system-dev libboost-thread-dev libsqlite3-dev libleveldb-dev leveldb-doc zlib1g-dev;
sudo apt-get -y install libboost-system-dev libboost-thread-dev libsqlite3-dev libtokyocabinet-dev zlib1g-dev;

#------------------------------------------------------------------------
mkdir -p $build_dir;mkdir $build_dir/percolator;mkdir $build_dir/converters;

######percolator########
#-----cmake-----
cd $build_dir/percolator;
echo -n "cmake percolator.....";
fakeroot -- cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr $src_dir/percolator;

#-----make------
echo -n "make percolator (this will take few minutes).....";
fakeroot -- make -j2;
fakeroot -- make -j2 package;

#######converters########
cd $build_dir/converters
#-----cmake-----
echo -n "cmake converters.....";
fakeroot -- cmake -DTARGET_ARCH=amd64 -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DSERIALIZE="TokyoCabinet" $src_dir/percolator/src/converters;

#-----make------
echo -n "make converters (this will take few minutes).....";

fakeroot -- make -j2;
fakeroot -- make -j2 package;

###########################
echo "Finished buildscript execution";
echo "in build directory ${build_dir}";
