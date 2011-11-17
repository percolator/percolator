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

###############################################################################
# compile and generate installers
rm -rf ${base}/percolatorBuild; mkdir ${base}/percolatorBuild; cd ${base}/percolatorBuild;
cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/percolator/cmake/windows_percolator/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-i686-windows;'$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $percolatorSource/percolator
make -j 8 win32installer


rm -rf ${base}/eludeBuild; mkdir ${base}/eludeBuild; cd ${base}/eludeBuild;
cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/percolator/src/elude/cmake/windows_elude/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $percolatorSource/percolator/src
make -j 8 win32installer

rm -rf ${base}/convertersBuild; mkdir ${base}/convertersBuild; cd ${base}/convertersBuild;
cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$percolatorSource/percolator/src/converters/cmake/windows_converter/Toolchain-mingw32.cmake '-DCMAKE_PREFIX_PATH='$base'/xsd-3.3.0-x86_64-linux-gnu' -DMINGW=ON -DSTATIC=ON -DGOOGLE_TEST=FALSE -DCMAKE_INSTALL_PREFIX=$base $percolatorSource/percolator/src
make -j 8 win32installer