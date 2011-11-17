#!/bin/sh
#TODO
# win: automatically turn on virtual machine
# askUser for interactive mode does not work: answering "n" once terminates the script

#run()
#{


set -e # as soon as a command fails terminate this script

if [ "$(id -u)" != "0" ]; then # check sudo rights
   echo "Error in $0 - this script must be run as root"
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
echo "Preparing to release version ${version_major}.${version_minor}..."

# set variables
base=/tmp/percolator_${version_major}_${version_minor}
percolatorSource=${base}/repository
packageDestination=${base}/packages

###############################################################################
# install libraries
installLibraries()
{
  echo "installing necessary libraries" 
  apt-get install libxerces-c-dev libboost-dev build-essential cmake libtokyocabinet-dev libsqlite3-dev zlib1g-dev libboost-filesystem-dev libboost-system-dev libbost-serialization-dev libleveldb-dev leveldb-doc
}

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
  mkdir ${percolatorSource}

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
  cd ${base}
  cs_pack="xsd-"${cs_version}".0-x86_64-linux-gnu" # xsd64
  download $cs_site/"x86_64" $cs_pack
  cs_pack="xsd-"${cs_version}".0-i686-linux-gnu" #xsd32
  download $cs_site/"i686" $cs_pack
  # clone repository
  cd ${percolatorSource}
  git clone git://github.com/percolator/percolator
}

###############################################################################
# change version number in the repository
changeVersion()
{
  echo "changing version in CMakeLists files" 
  # change version number for percolator
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MAJOR "[0-9]*}{CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'}g' ${percolatorSource}/percolator/CMakeLists.txt
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MINOR "[0-9]*}{CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'}g' ${percolatorSource}/percolator/CMakeLists.txt
  # change version number for the converters
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MAJOR "[0-9]*}{CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'}g' ${percolatorSource}/percolator/src/converters/CMakeLists.txt
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MINOR "[0-9]*}{CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'}g' ${percolatorSource}/percolator/src/converters/CMakeLists.txt
  # change version number for elude
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MAJOR "[0-9]*}{CPACK_PACKAGE_VERSION_MAJOR "'${version_major}'}g' ${percolatorSource}/percolator/src/elude/CMakeLists.txt
  perl -pi -e 's{CPACK_PACKAGE_VERSION_MINOR "[0-9]*}{CPACK_PACKAGE_VERSION_MINOR "'${version_minor}'}g' ${percolatorSource}/percolator/src/elude/CMakeLists.txt
}

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
percolator64()
{
  object="PERCOLATOR64bit"
  echo ${object}
  rm -rf ${base}/percolatorBuild; mkdir ${base}/percolatorBuild; cd ${base}/percolatorBuild;
  cmake -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/percolator
  make -j 8
  make package
  make package_source
  rm -rf ${packageDestination}/percolator64; mkdir ${packageDestination}/percolator64; # harvest packages
  cp *.deb ${packageDestination}/percolator64
  cp *.rpm ${packageDestination}/percolator64
  cp *.tar.gz ${packageDestination}
}

###############################################################################
converters64()
{
  object="CONVERTERS64bit"
  echo ${object}
  rm -rf ${base}/percolator-convertersBuild; mkdir ${base}/percolator-convertersBuild; cd ${base}/percolator-convertersBuild;
  cmake -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/percolator/src/converters
  make -j 8
  make package
  rm -rf ${packageDestination}/converters64; mkdir ${packageDestination}/converters64; # harvest packages
  cp *.deb ${packageDestination}/converters64
  cp *.rpm ${packageDestination}/converters64
}

###############################################################################
elude64()
{
  object="ELUDE64bit"
  echo ${object}
  rm -rf ${base}/eludeBuild; mkdir ${base}/eludeBuild; cd ${base}/eludeBuild;
  cmake -DCMAKE_BUILD_TYPE=Release '-DCMAKE_PREFIX_PATH='${base}'/xsd-3.3.0-x86_64-linux-gnu/' -DCMAKE_INSTALL_PREFIX='/usr' ${percolatorSource}/percolator/src/elude
  make -j 8
  make package
  rm -rf ${packageDestination}/elude64; mkdir ${packageDestination}/elude64; # harvest packages
  cp *.deb ${packageDestination}/elude64
  cp *.rpm ${packageDestination}/elude64
}

###############################################################################
windows()
{
  object="WINDOWS"
  echo ${object}
  ssh -l root -p 2222 localhost 'bash -s' < ${percolatorSource}/percolator/install/release_windows.sh ${version_major} ${version_minor} && connection="open" || connection="closed"
  if [ "$connection" = "open" ]; then
    # harvest packages
    rm -rf ${packageDestination}/win; mkdir ${packageDestination}/win;
    scp -P 2222 root@localhost:/tmp/percolator_2_03/percolatorBuild/*-win32.exe ${packageDestination}/win
    scp -P 2222 root@localhost:/tmp/percolator_2_03/eludeBuild/*-win32.exe ${packageDestination}/win
    scp -P 2222 root@localhost:/tmp/percolator_2_03/convertersBuild/*-win32.exe ${packageDestination}/win
  else
    echo "The remote virtual machine is offline. Skipping windows executables."
  fi
}


###############################################################################
installLibraries
clean
changeVersion
percolator64
converters64
elude64
windows

echo ""
echo "Release process for version ${version_major}.${version_minor} 64 bit completed successfully!"
echo "The newly generated packages can be found in ${packageDestination}"

#} # end of run()
#if [ -f /tmp/percolator_release_error.log ]; then
#  rm /tmp/percolator_release_error.log
#fi
#run $1 $2 2> /tmp/percolator_release_error.log
