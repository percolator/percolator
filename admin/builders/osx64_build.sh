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
    other_packages="cmake tokyocabinet bzip2 libiconv zlib gtest"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    boost_install_options="boost"
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib libomp googletest"
elif [[ -f /usr/local/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager="brew"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    boost_install_options="boost"
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib libomp googletest"

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
	curl -k -O ${mac_os_xerces_url}
	tar xzf ${mac_os_xerces}.tar.gz
	cd ${mac_os_xerces}/
	./configure CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64" --disable-dynamic --enable-transcoder-iconv --disable-network --disable-threads
	make -j 2
	sudo make install
fi

cd ${build_dir}

# XSD installation
if [ ! -d ${mac_os_xsd} ]; then
#  if [ $package_manager == "sudo port" ]; then
#     export XSDDIR=/usr/local/Cellar/xsd/4.0.0_1/
#  fi
  curl -OL ${mac_os_xsd_url}
  tar -xjf ${mac_os_xsd}.tar.bz2
  cd ${mac_os_xsd}
  
  # https://www.codesynthesis.com/pipermail/xsde-users/2022-August/000916.html
  mv libxsd-frontend/version libxsd-frontend/version.txt
  mv libcutl/version libcutl/version.txt
  mv xsd/version xsd/version.txt

  echo '#include <iostream>' > tmp_file
  cat libxsd-frontend/xsd-frontend/semantic-graph/elements.cxx >>  tmp_file
  mv tmp_file libxsd-frontend/xsd-frontend/semantic-graph/elements.cxx
  
  make CPPFLAGS=-I../${mac_os_xerces}/src LDFLAGS=-L../${mac_os_xerces}/src/.libs
  ./xsd/xsd/xsd --version
  # Move Binary, to the right include files
  mv xsd/xsd/xsd xsd/libxsd/xsd/
  rm -fr xsd/xsd
  cd ..
  export XSDDIR=${build_dir}/${mac_os_xsd}/xsd/libxsd
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

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DXML_SUPPORT=ON -DGOOGLE_TEST=1 -DCMAKE_PREFIX_PATH="${build_dir}/${mac_os_xerces}/;${build_dir}/${mac_os_xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/"  ${src_dir}/percolator
make -j 2
make -j 2 package

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters

cmake -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_PREFIX_PATH="${build_dir}/${mac_os_xerces}/;${build_dir}/${mac_os_xsd}/;/opt/local/;/usr/;/usr/local/;~/;/Library/Developer/CommandLineTools/usr/" -DSERIALIZE="TokyoCabinet" ${src_dir}/percolator/src/converters
make -j 2
make -j 2 package
#--------------------------------------------

echo "build directory was : ${build_dir}";

cp -v ${build_dir}/{percolator-noxml,percolator,converters}/*.pkg ${release_dir};
