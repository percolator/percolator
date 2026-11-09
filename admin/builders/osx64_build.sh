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
    other_packages="cmake bzip2 libiconv zlib"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    other_packages="pbzip2 lzlib llvm libomp"
elif [[ -f /usr/local/bin/brew || -f /opt/homebrew/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager="brew"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    other_packages="pbzip2 lzlib llvm libomp"
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

if command -v cmake >/dev/null 2>&1; then
  echo "CMake found, stay with installation"
else
  other_packages="$other_packages cmake"
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

if [[ -f /usr/local/bin/brew || -f /opt/homebrew/bin/brew ]]
  then
  export HOMEBREW_PREFIX="$(brew --prefix)"
  export OMP_PREFIX="$HOMEBREW_PREFIX/opt/libomp"

  # Make omp.h and libomp visible
  export CXXFLAGS="$CXXFLAGS -I$OMP_PREFIX/include"
  export LDFLAGS="$LDFLAGS -L$OMP_PREFIX/lib"
fi

cd ${src_dir}

# read all urls and file names from a centralized kb file
source ./percolator/admin/builders/_urls_and_file_names_.sh
mkdir -p ${build_dir}
cd ${build_dir}

function build_component() {
    local name="$1"
    local src="$2"
    local options="$3"

    mkdir -p "${build_dir}/${name}"
    pushd "${build_dir}/${name}" > /dev/null

    cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr/local \
      -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
      -DCMAKE_EXE_LINKER_FLAGS="$LDFLAGS" \
      $options "$src"
    make -j$(sysctl -n hw.ncpu)
    make package
    popd > /dev/null
}

mkdir -p ${release_dir}

build_component "percolator" "${src_dir}/percolator" "-DGOOGLE_TEST=1"

echo "build directory was : ${build_dir}";

cp -v ${build_dir}/percolator/*.pkg ${release_dir};
