#!/bin/bash
# Requirements are:
# Command line tools
# MacPorts or homebrew as package manager
# CMake
# PackageMaker (https://developer.apple.com/downloads ::search::
#----------------------------------------


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

if [[ ! -d /Applications/XCode.app ]]
  then
    echo "Apple developer tools are required (Search for XCode in the App Store)"
    exit 1
fi

package_manager_installed=true
if [[ -d /opt/local/var/macports ]]
  then
    echo "[ Package manager ] : MacPorts "
    package_manager="sudo port"
    boost_install_options="boost -no_static -no_single"
    other_packages="tokyocabinet bzip2 libiconv zlib pthread"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    boost_install_options="boost"
    other_packages="tokyo-cabinet lbzip2 pbzip2 lzlib"
else
    package_manager_installed=false
fi

if [ "$package_manager_installed" == false ]
  then
  echo "Error: no suitable package manager installed"
  echo " Get homebrew or macports:"
  echo "  Homebrew: http://brew.sh/ "
  echo "  MacPorts: http://www.macports.org/install.php"
  exit 1
fi


# PackageMaker is also required. It is a part of Auxiliary Tools for Xcode packet that Apple provides for Developers.
# You can find it here: https://developer.apple.com/downloads/index.action?name=packagemaker#
# TODO: check if PackageMaker is installed

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d -t build)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    if [[ ! -f /usr/bin/git ]]; then
      $package_manager install git;
    fi
    src_dir="$(mktemp -d -t src)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
	src_dir="${src_dir}/percolator"
  else
    # Might not work if we have symlinks in the way
    src_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../" && pwd )
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for the user"
whoami

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for user" `whoami`
# exit 1
$package_manager install $other_packages
$package_manager install $boost_install_options

#----------------------------------------
xsd=xsd-3.3.0-i686-macosx
cd ${build_dir}
curl -O http://www.codesynthesis.com/download/xsd/3.3/macosx/i686/${xsd}.tar.bz2
tar -xjf ${xsd}.tar.bz2
sed -i -e 's/setg/this->setg/g' ${xsd}/libxsd/xsd/cxx/zc-istream.txx
sed -i -e 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/parsing.txx
sed -i -e 's/ push_back/ this->push_back/g' ${xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx

extr=${xsd}/libxsd/xsd/cxx/tree/xdr-stream-extraction.hxx
inse=${xsd}/libxsd/xsd/cxx/tree/xdr-stream-insertion.hxx
sed -i -e 's/ uint8_t/ unsigned char/g' ${extr} ${inse}
sed -i -e 's/ int8_t/ char/g' ${extr} ${inse}
sed -i -e 's/xdr_int8_t/xdr_char/g' ${extr} ${inse}
sed -i -e 's/xdr_uint8_t/xdr_u_char/g' ${extr} ${inse}
#------------------------------------------
xer=xerces-c-3.1.1
mkdir -p ${build_dir}
cd ${build_dir}
if [[ -d /usr/local/include/xercesc ]]
	then
	echo "Xerces is already installed."
else
	curl -O http://apache.mirrors.spacedump.net/xerces/c/3/sources/${xer}.tar.gz
	tar xzf ${xer}.tar.gz
	cd ${xer}/
	./configure CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64" --disable-dynamic --enable-transcoder-iconv --disable-network --disable-threads
	make -j 2
	sudo make install
fi

#-------------------------------------------

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml

cmake -DCXX="/usr/bin/gcc" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/ -DXML_SUPPORT=OFF -DCMAKE_PREFIX_PATH="${build_dir}/${xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  ${src_dir}/
make -j 2 "VERBOSE=1"
sudo make -j 2 package


mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

cmake -DCXX="/usr/bin/gcc" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/ -DXML_SUPPORT=ON -DCMAKE_PREFIX_PATH="${build_dir}/${xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  ${src_dir}/
make -j 2 "VERBOSE=1"
sudo make -j 2 package

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters

cmake -V  -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/ -DCMAKE_PREFIX_PATH="${build_dir}/${xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  -DSERIALIZE="TokyoCabinet" ${src_dir}/src/converters
make -j 2 "VERBOSE=1"
sudo make -j 2 package
#--------------------------------------------
mkdir -p ${release_dir}
cp -v ${build_dir}/percolator-noxml/*.dmg ${release_dir}
cp -v ${build_dir}/percolator/*.dmg ${release_dir}
cp -v ${build_dir}/converters/*.dmg ${release_dir}
