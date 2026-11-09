#!/bin/bash

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

if [[ -z ${build_dir} ]]; then
  build_dir="${HOME}/build";
#  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir src_XXXX)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for the user"


sudo yum install -y cmake wget mingw-w64-tools mingw32-gcc-c++ mingw32-filesystem mingw-binutils-generic mingw32-nsis
sudo yum install -y mingw32-sqlite mingw32-zlib mingw32-curl mingw32-pthreads

cd ${src_dir}

mkdir -p ${build_dir}
cd ${build_dir}

# download, compile and link percolator

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator

mingw32-cmake -DCMAKE_BUILD_TYPE=Release ${src_dir}/percolator
make -j4 package

cp -v per*.exe ${release_dir}

echo "build directory is : ${build_dir}";
