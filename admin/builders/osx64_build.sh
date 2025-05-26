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

# get current architecture
ARCH=$(uname -m)

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
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib llvm libomp"
elif [[ -f /usr/local/bin/brew || -f /opt/homebrew/bin/brew ]]  
  then
    echo "[ Package manager ] : Homebrew "
    package_manager="brew"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    boost_install_options="boost"
    other_packages="cmake tokyo-cabinet lbzip2 pbzip2 lzlib llvm libomp"

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

$package_manager install xerces-c xsd
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
export SYSROOT_FLAGS="--sysroot=$SDKROOT -isystem$SDKROOT/usr/include -isystem$SDKROOT/System/Library/Frameworks"
export HOMEBREW_LLVM_INCLUDE="/opt/homebrew/opt/llvm/include/c++/v1"
export CXXFLAGS="-isystem $HOMEBREW_LLVM_INCLUDE $SYSROOT_FLAGS -I/opt/homebrew/opt/tokyo-cabinet/include"
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib/c++ -L/opt/homebrew/opt/tokyo-cabinet/lib -stdlib=libc++ -framework CoreFoundation -framework CoreServices -lcurl"
prefix_path="/opt/homebrew/opt/xerces-c;/opt/homebrew/opt/xsd"


function build_component() {
    local name="$1"
    local src="$2"
    local options="$3"

    mkdir -p "${build_dir}/${name}"
    pushd "${build_dir}/${name}" > /dev/null

    cmake \
      -DCMAKE_C_COMPILER="$CC" \
      -DCMAKE_CXX_COMPILER="$CXX" \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr/local \
      -DCMAKE_PREFIX_PATH="$prefix_path" \
      -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
      -DCMAKE_EXE_LINKER_FLAGS="$LDFLAGS" \
      $options "$src"
    make -j$(sysctl -n hw.ncpu)
    make package
    popd > /dev/null
}

mkdir -p ${release_dir}

build_component "percolator-noxml" "${src_dir}/percolator" "-DXML_SUPPORT=OFF"
build_component "percolator" "${src_dir}/percolator" "-DXML_SUPPORT=ON -DGOOGLE_TEST=1"
build_component "converters" "${src_dir}/percolator/src/converters" "-DSERIALIZE=TokyoCabinet"

echo "build directory was : ${build_dir}";

cp -v ${build_dir}/{percolator-noxml,percolator,converters}/*.pkg ${release_dir};
