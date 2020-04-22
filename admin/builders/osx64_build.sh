#!/bin/bash
# Requirements are:
# XCode
# Command line tools (check if installed with "xcode-select -p", otherwise install with "xcode-select --install")
# MacPorts or homebrew as package manager
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

#if [[ ! -d /Applications/PackageMaker.app ]]
#  then
#    echo "Apple developer PackageManager is required and expected in the "
#    echo "/Applications folder. If you have moved it elsewhere, please change this script"
#    echo ""
#    echo "It is part of the Auxiliary tools for XCode - Late July 2012"
#    echo "Yes, 2012! since then Apple moved to the app store and requires"
#    echo "packages and dmgs to be build differently. "
#    echo "However, the old packagemaker still works with 10.11"
#    echo
#    echo "You can find it here: "
#    echo "http://adcdownload.apple.com/Developer_Tools/auxiliary_tools_for_xcode__late_july_2012/xcode44auxtools6938114a.dmg"
#    echo ""
#    exit 1
#fi

package_manager_installed=true
if [[ -d /opt/local/var/macports ]]
  then
    echo "[ Package manager ] : MacPorts "
    package_manager="sudo port"
    boost_install_options="boost -no_static"
    other_packages="cmake tokyocabinet bzip2 libiconv zlib"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    boost_install_options="boost"
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib xerces-c xsd"
elif [[ -f /usr/local/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager="brew"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    boost_install_options="boost"
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib xerces-c xsd"

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
    src_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../" && pwd )
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for user" `whoami`
$package_manager install $other_packages
$package_manager install $boost_install_options

cd ${src_dir}

# read all urls and file names from a centralized kb file
source ./percolator/admin/builders/_urls_and_file_names_.sh
mkdir -p ${build_dir}
cd ${build_dir}

# XercesC installation
if [[ -d /usr/local/include/xercesc ]] # this implies homebrew installation ...
	then
	echo "Xerces is already installed."
else
	curl -O ${mac_os_xerces_url}
	tar xzf ${mac_os_xerces}.tar.gz
	cd ${mac_os_xerces}/
	./configure CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64" --disable-dynamic --enable-transcoder-iconv --disable-network --disable-threads
	make -j 2
	sudo make install
fi

cd ${build_dir}

# XSD installation
if [ ! -d ${mac_os_xsd} ] && [ $package_manager == "sudo port" ]; then
  curl -OL ${mac_os_xsd_url}
  tar -xjf ${mac_os_xsd}.tar.bz2
  sed -i -e 's/setg/this->setg/g' ${mac_os_xsd}/libxsd/xsd/cxx/zc-istream.txx
  sed -i -e 's/ push_back/ this->push_back/g' ${mac_os_xsd}/libxsd/xsd/cxx/tree/parsing.txx
  sed -i -e 's/ push_back/ this->push_back/g' ${mac_os_xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx

  extr=${mac_os_xsd}/libxsd/xsd/cxx/tree/xdr-stream-extraction.hxx
  inse=${mac_os_xsd}/libxsd/xsd/cxx/tree/xdr-stream-insertion.hxx
  sed -i -e 's/ uint8_t/ unsigned char/g' ${extr} ${inse}
  sed -i -e 's/ int8_t/ char/g' ${extr} ${inse}
  sed -i -e 's/xdr_int8_t/xdr_char/g' ${extr} ${inse}
  sed -i -e 's/xdr_uint8_t/xdr_u_char/g' ${extr} ${inse}
else
  echo "XSD is already installed."
fi

#-------------------------------------------

mkdir -p ${release_dir}

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DXML_SUPPORT=OFF -DCMAKE_PREFIX_PATH="/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  ${src_dir}/percolator
make -j 2
make -j 2 package
cp -v per*.pkg ${release_dir}

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DXML_SUPPORT=ON -DCMAKE_PREFIX_PATH="${build_dir}/${mac_os_xerces}/;${build_dir}/${mac_os_xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  ${src_dir}/percolator
make -j 2
make -j 2 package
cp -v per*.pkg ${release_dir}

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_PREFIX_PATH="${build_dir}/${mac_os_xerces}/;${build_dir}/${mac_os_xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/" -DSERIALIZE="TokyoCabinet" ${src_dir}/percolator/src/converters
make -j 2
make -j 2 package
cp -v per*.pkg ${release_dir}

mkdir -p ${build_dir}/elude
cd ${build_dir}/elude

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_PREFIX_PATH="${build_dir}/${mac_os_xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/" ${src_dir}/percolator/src/elude_tool
make -j 2
make -j 2 package
cp -v elude*.pkg ${release_dir}
#--------------------------------------------

echo "build directory was : ${build_dir}";
