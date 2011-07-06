#!/bin/sh
#TODO
# find more robust way of setting version numbers: current assumption is that version is set to 0.00
# revert changes to reporitory (permissions to CMakeLists and new version numbers)
# win: find a way to execute script throught ssh in batch mode (remember it needs to be run from the root of the repo)
# win: automatically turn on virtual machine
# duplicate instead of redirecting stderr (run() + last few lines)
# askUser for interactive mode does not work: answering "n" once terminates the script
# 32 bit versions: CMAKE_CXX_FLAGS="-m32" broken

#run()
#{


set -e # as soon as a command fails terminate this script

if [ "$(id -u)" != "0" ]; then
   echo "Error in $0 - this script must be run as root"
   exit 1
fi
if ! [ -f CMakeLists.txt ]; then # check script was invoked from right location
  echo "Error in $0 - please invoke the script from the root of the repository"
  exit 1
fi
# get command line parameters
if [ $# -ne 2 ]; then # check number of parameters
  echo "Error in $0 - Invalid Argument Count"
  echo "Syntax: $0 [version_major] [version_minor]"
  exit 1
fi
version_major=$1
version_minor=$2
echo "executing release scritp for version ${version_major}.${version_minor}"

# set variables
base=/tmp/percolator_${version_major}_${version_minor}
percolatorSource=`pwd`
packageDestination=${base}/packages
cd ${base}

# install libraries
echo "installing necessary libraries" 
apt-get install libxerces-c-dev libboost-dev build-essential cmake libtokyocabinet-dev libsqlite3-dev

###############################################################################
# clean traces of past runs of this script
clean()
{
  echo "cleaning working directory" ${base}
  if [ -d "${base}" ]; then
    rm -r ${base}
  fi
  mkdir ${base}
  mkdir ${packageDestination}
  cd ${base}

  download(){
    site=$1
    pack=$2
    wget -nc "${site}/${pack}.tar.bz2"
    tar xjf "${pack}.tar.bz2"
    rm "${pack}.tar.bz2"
  }
  # install xsd
  cs_version="3.3"
  cs_site="http://codesynthesis.com/download/xsd/"${cs_version}"/linux-gnu"
  echo "downloading codesynthesis version"  ${cs_version}
  cs_pack="xsd-"${cs_version}".0-x86_64-linux-gnu" # xsd64
  download $cs_site/"x86_64" $cs_pack
  cs_pack="xsd-"${cs_version}".0-i686-linux-gnu" #xsd32
  download $cs_site/"i686" $cs_pack
}
clean

###############################################################################
# change version number in the repository
changeVersion()
{
  # change version number for percolator
  sed 's/CPACK_PACKAGE_VERSION_MAJOR "0/CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'/' <${percolatorSource}/CMakeLists.txt >${percolatorSource}/CMakeLists_tmp.txt
  rm ${percolatorSource}/CMakeLists.txt
  sed 's/CPACK_PACKAGE_VERSION_MINOR "00/CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'/' <${percolatorSource}/CMakeLists_tmp.txt >${percolatorSource}/CMakeLists.txt
  rm ${percolatorSource}/CMakeLists_tmp.txt
  chown tomasoni ${percolatorSource}/CMakeLists.txt # TODO this is a hack
  chgrp users ${percolatorSource}/CMakeLists.txt
  # change version number for the converters
  sed 's/CPACK_PACKAGE_VERSION_MAJOR "0/CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'/' <${percolatorSource}/src/converters/CMakeLists.txt >${percolatorSource}/src/converters/CMakeLists_tmp.txt
  rm ${percolatorSource}/src/converters/CMakeLists.txt
  sed 's/CPACK_PACKAGE_VERSION_MINOR "00/CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'/' <${percolatorSource}/src/converters/CMakeLists_tmp.txt >${percolatorSource}/src/converters/CMakeLists.txt
  rm ${percolatorSource}/src/converters/CMakeLists_tmp.txt
  chown tomasoni ${percolatorSource}/src/converters/CMakeLists.txt # TODO this is a hack
  chgrp users ${percolatorSource}/src/converters/CMakeLists.txt
  # change version number for elude
  sed 's/CPACK_PACKAGE_VERSION_MAJOR "0/CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'/' <${percolatorSource}/src/elude/CMakeLists.txt >${percolatorSource}/src/elude/CMakeLists_tmp.txt
  rm ${percolatorSource}/src/elude/CMakeLists.txt
  sed 's/CPACK_PACKAGE_VERSION_MINOR "00/CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'/' <${percolatorSource}/src/elude/CMakeLists_tmp.txt >${percolatorSource}/src/elude/CMakeLists.txt
  rm ${percolatorSource}/src/elude/CMakeLists_tmp.txt
  chown tomasoni ${percolatorSource}/src/elude/CMakeLists.txt # TODO this is a hack
  chgrp users ${percolatorSource}/src/elude/CMakeLists.txt
}
changeVersion

#askUser(){ #TODO interactive mode, can't get it to work!
#  echo "Generate packages for" $1"? [y/n]"; 
#  read answer
#  if [ "${answer}" = "y" ]; then
#    return "0"
#  else
#    return "1"
#  fi
#}
#askUser "percolator64"
#if [ "$?" = "0" ]; then
#  percolator64
#fi


###############################################################################
object="PERCOLATOR64bit"
percolator64()
{
  echo ${object}
  rm -rf ${base}/percolatorBuild; mkdir ${base}/percolatorBuild; cd ${base}/percolatorBuild;
  cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}
  make -j 8
  make package
  make package_source
  rm -rf ${packageDestination}/percolator64; mkdir ${packageDestination}/percolator64; # harvest packages
  cp *.deb ${packageDestination}/percolator64
  cp *.rpm ${packageDestination}/percolator64
  cp *.tar.gz ${packageDestination}
}
percolator64

###############################################################################
object="PERCOLATOR32bit"
percolator32()
{
  echo ${object}
  rm -rf ${base}/percolatorBuild; mkdir ${base}/percolatorBuild; cd ${base}/percolatorBuild;
  cmake -DTARGET_ARCH=i386 -DCMAKE_CXX_FLAGS="-m32" -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-i686-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}
  make -j 8
  make package
  rm -rf ${packageDestination}/percolator32; mkdir ${packageDestination}/percolator32; # harvest packages
  cp *.deb ${packageDestination}/percolator32
  cp *.rpm ${packageDestination}/percolator32
}
#percolator32

###############################################################################
object="CONVERTERS64bit"
converters64()
{
  echo ${object}
  rm -rf ${base}/percolator-convertersBuild; mkdir ${base}/percolator-convertersBuild; cd ${base}/percolator-convertersBuild;
  cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/src/converters
  make -j 8
  make package
  rm -rf ${packageDestination}/converters64; mkdir ${packageDestination}/converters64; # harvest packages
  cp *.deb ${packageDestination}/converters64
  cp *.rpm ${packageDestination}/converters64
}
converters64

###############################################################################
object="CONVERTERS32bit"
converters32()
{
  echo ${object}
  rm -rf ${base}/percolator-convertersBuild; mkdir ${base}/percolator-convertersBuild; cd ${base}/percolator-convertersBuild;
  cmake -DTARGET_ARCH=i386 -DCMAKE_CXX_FLAGS="-m32" -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-i686-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/src/converters
  make -j 8
  make package
  rm -rf ${packageDestination}/converters32; mkdir ${packageDestination}/converters32; # harvest packages
  cp *.deb ${packageDestination}/converters32
  cp *.rpm ${packageDestination}/converters32
}
#converters32

###############################################################################
object="ELUDE64bit"
elude64()
{
  echo ${object}
  rm -rf ${base}/eludeBuild; mkdir ${base}/eludeBuild; cd ${base}/eludeBuild;
  cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/src/elude
  make -j 8
  make package
  rm -rf ${packageDestination}/elude64; mkdir ${packageDestination}/elude64; # harvest packages
  cp *.deb ${packageDestination}/elude64
  cp *.rpm ${packageDestination}/elude64
}
elude64

###############################################################################
object="ELUDE32bit"
elude32()
{
  echo ${object}
  rm -rf ${base}/eludeBuild; mkdir ${base}/eludeBuild; cd ${base}/eludeBuild;
  cmake -DTARGET_ARCH=i386 -DCMAKE_CXX_FLAGS="-m32" -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-i686-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/src/elude
  make -j 8
  make package
  rm -rf ${packageDestination}/elude32; mkdir ${packageDestination}/elude32; # harvest packages
  cp *.deb ${packageDestination}/elude32
  cp *.rpm ${packageDestination}/elude32
}
#elude32

###############################################################################
#object="WINDOWS"
#windows()
#{
#  ssh -l root -p 2222 localhost #TODO this opens an interactive session! try to do it in batch mode
#  cd /mnt/VirtualBoxShare/percolator
#  ./install/windows_percolator.sh
#  cd src/elude
#  ./install/windows_elude.sh
#  logout
#  rm -rf ${packageDestination}/win; mkdir ${packageDestination}/win; # harvest packages
#  cp /scratch/VirtualBoxShare/percolatorBuild/*.exe ${packageDestination}/win
#  cp /scratch/VirtualBoxShare/eludeBuild/*.exe ${packageDestination}/win
#}
#windows


#} # end of run()
#if [ -f /tmp/percolator_release_error.log ]; then
#  rm /tmp/percolator_release_error.log
#fi
#run $1 $2 2> /tmp/percolator_release_error.log
