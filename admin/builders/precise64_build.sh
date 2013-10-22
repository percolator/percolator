#!/bin/bash

# managing input arguments

if [ $# -eq 2 ];
  then src_dir=$1;build_dir=$2;
elif [ $# -eq 1 ];
  then sudo apt-get install git;
  tmp_dir="$(mktemp -d --tmpdir precise_tmp_XXXX)";
  mkdir "$tmp_dir"/src;mkdir "$tmp_dir"/build;
  src_dir=""$tmp_dir"/src";build_dir=""$tmp_dir"/build";
  git clone --branch "$1" https://github.com/percolator/percolator.git "$src_dir"/percolator;
else 
  echo "Please add either one argument as branch name or two arguments for your source directory containing percolator/ and build directory";
  exit 1;
fi;

echo "Building the Percolator packages with src=${src} and build=${build}"


#------------------------------------------------------------------------
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

sudo apt-get update;
sudo apt-get -y install g++ make cmake rpm;
sudo apt-get -y install xsdcxx libxerces-c-dev libboost-dev libboost-filesystem-dev;
sudo apt-get -y install libboost-system-dev libboost-thread-dev libsqlite3-dev libleveldb-dev leveldb-doc zlib1g-dev;

#------------------------------------------------------------------------
mkdir -p $build_dir;mkdir $build_dir/percolator;mkdir $build_dir/converters;

######percolator########
#-----cmake-----
cd $build_dir/percolator;
echo -n "cmake percolator.....";
if cmake -DTARGET_ARCH=amd64 -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" $src_dir/percolator;
then echo "Done";
else echo "Cmake was unsuccessful!";fi;
#-----make------
echo -n "make percolator (this will take few minutes).....";
if make -j2 package;
then echo "Done";
else echo "make was unsuccessful!";fi;
#######converters########
cd $build_dir/converters
#-----cmake-----
echo -n "cmake converters.....";
if cmake -DTARGET_ARCH=amd64 -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" $src_dir/percolator/src/converters;
then echo "Done";
else echo "Cmake was unsuccessful!";fi;
#-----make------
echo -n "make converters (this will take few minutes).....";
if make -j2 package;
then echo "Done";
else echo "make was unsuccessful!";fi;
###########################
echo "build directory was : "$build_dir"";
