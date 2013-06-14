#!/bin/bash

branch="branch-2-05";
src_dir="percolator";

function pkg_mng()
{
	echo -n "installing the package $1 ...";
	(sudo apt-get -y install $1) > /dev/null;
	if [ $? -eq 0 ] 
	then echo "$1 is installed."; return 0;
	else echo "$1 could not be installed."; return 1;
	fi;
}
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

sudo apt-get update;
pkg_mng "g++";
pkg_mng "make";
pkg_mng "cmake";
pkg_mng "rpm";
pkg_mng "git";
pkg_mng "xsdcxx";
pkg_mng "libxerces-c-dev";
pkg_mng "libboost-dev";
pkg_mng "libboost-filesystem-dev";
pkg_mng "libboost-system-dev";
pkg_mng "libboost-thread-dev";
pkg_mng "libsqlite3-dev";
pkg_mng "libleveldb-dev";
pkg_mng "leveldb-doc";
pkg_mng "zlib1g-dev";

#------------------------------------------------------------------------
cd /vagrant; 
cur_dir=$PWD; #if there is no /vagrant
#-----clone-----
rm -rf $src_dir;
echo -n "Cloning Percolator branch:$branch.....";
 git clone https://github.com/percolator/percolator.git --branch $branch;
######percolator########
#-----cmake-----
rm -rf build;mkdir build;mkdir build/percolator;mkdir build/converters;
cd $cur_dir/build/percolator;
echo -n "cmake percolator.....";
if (cmake $cur_dir/$src_dir) > /dev/null ;
then echo "Done";
else echo "Cmake was unsuccessful!";return 1; fi;
#-----make------
echo -n "make percolator (this will take few minutes).....";
if (make -j2 package) > /dev/null;
then echo "Done";
else echo "make was unsuccessful!";return 1; fi;
#######converters########
cd $cur_dir/build/converters
#-----cmake-----
echo -n "cmake converters.....";
if (cmake $cur_dir/$src_dir/src/converters) > /dev/null ;
then echo "Done";
else echo "Cmake was unsuccessful!";return 1; fi;
#-----make------
echo -n "make converters (this will take few minutes).....";
if (make -j2 package) > /dev/null;
then echo "Done";
else echo "make was unsuccessful!";return 1; fi;
###########################
#-----release-----
cd $cur_dir;mkdir release;
cp $cur_dir/build/percolator/per*.rpm release/
cp $cur_dir/build/percolator/per*.deb release/
cp $cur_dir/build/converters/per*.rpm release/
cp $cur_dir/build/converters/per*.deb release/
