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
  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir build_XXXX)";
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} and build=${build_dir} for the user"
whoami;

sudo dnf install -y gcc gcc-c++ make cmake wget rpm-build
sudo dnf install -y tokyocabinet-devel boost-static boost-devel sqlite-devel zlib-devel bzip2-devel xerces-c-devel xsd

cd ${src_dir}

# download, compile and link percolator

mkdir -p ${build_dir}/percolator-noxml
cd ${build_dir}/percolator-noxml
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF ${src_dir}/percolator
make -j 4;
make -j 4 package;

mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=ON ${src_dir}/percolator
make -j 4;
make -j 4 package;

mkdir -p ${build_dir}/converters
cd ${build_dir}/converters
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DSERIALIZE="TokyoCabinet" ${src_dir}/percolator/src/converters
make -j 4;
make -j 4 package;

mkdir -p ${build_dir}/elude
cd ${build_dir}/elude
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr ${src_dir}/percolator/src/elude_tool
make -j 4;
make -j 4 package;

echo "build directory was : ${build_dir}";

cp -v ${build_dir}/{percolator-noxml,percolator,converters,elude}/*.rpm ${release_dir};

