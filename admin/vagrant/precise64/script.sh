#!/bin/bash

branch="branch-2-05";

src_dir="percolator";
function pkg_check()
{
	echo -n "$1....."
	(dpkg -s $1) &> /dev/null;
	c=$?;
	if [ $c -eq 0 ]
	then
		echo  "OK";
		return 0;
	else 
		echo  "NOT installed";
		return 1; 
	fi
}

function pkg_install()
{
	echo -n "installing the package $1 ...";
	(sudo apt-get -y install $1) > /dev/null;
	if [ $? -eq 0 ] 
	then echo "$1 is installed successfully."; return 0;
	else echo "$1 could not be installed."; return 1;
	fi;
}

function pkg_mng()
{
	pkg_check $1;
	if [ $? -eq 1 ]; then pkg_install $1; fi;
	return;
}
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

pkg_mng "g++";
pkg_mng "make";
pkg_mng "cmake";
pkg_mng "rpm"; #needed for make -j2
pkg_mng "git";
pkg_mng "xsdcxx";
pkg_mng "libxerces-c-dev";
pkg_mng "libboost-dev";
pkg_mng "libboost-filesystem-dev";
pkg_mng "libboost-system-dev";
pkg_mng "libboost-thread-dev";
#------------------------------------------------------------------------
cd /vagrant; 
#-----clone-----
rm -rf $src_dir;
echo -n "Cloning Percolator branch:$branch.....";
 git clone git://github.com/percolator/percolator --branch $branch;
#-----cmake-----
rm -rf build;mkdir build;cd build;
echo -n "cmake.....";
if (cmake ../$src_dir) > /dev/null ;
then echo "Done";
else echo "Cmake was unsuccessful!";return 1; fi;
#-----make------
echo -n "make (this will take few minutes).....";
if (make -j2 package) > /dev/null;
then echo "Done";
else echo "make was unsuccessful!";return 1; fi;
#------------------------------------------------------------------------


