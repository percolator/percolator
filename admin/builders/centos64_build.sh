#!/bin/bash

[ "$(id -u)" -eq 0 ] && SUDO="" || SUDO="sudo"

# managing input arguments
while getopts "s:b:r:t:" OPTION; do
  case $OPTION in
    s) src_dir=${OPTARG};;
    t) branch=${OPTARG};;
    r) release_dir=${OPTARG};;
    b) build_dir=${OPTARG};;
    \?) echo "Invalid option: -${OPTARG}" >&2;;
  esac
done

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d --tmpdir build_XXXX)"
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    ${SUDO} dnf -y install git
    src_dir="$(mktemp -d --tmpdir build_XXXX)"
    git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator"
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} and build=${build_dir} for the user"
whoami

${SUDO} dnf clean all
${SUDO} dnf -y distro-sync
${SUDO} dnf -y install gcc gcc-c++ wget rpm-build cmake git
${SUDO} dnf -y install zlib-devel bzip2-devel
${SUDO} dnf -y install epel-release
${SUDO} dnf -y install gtest-devel

cd ${src_dir}

# download, compile and link percolator

echo "Installing percolator"

echo "cmake percolator ....."
mkdir -p ${build_dir}/percolator
cd ${build_dir}/percolator
(set -x; \
cmake -DTARGET_ARCH=x86_64 -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=1 -DCMAKE_INSTALL_PREFIX=/usr ${src_dir}/percolator; \
)
make -j 4
make test
make -j 4 package

echo "build directory was : ${build_dir}"

cp -v ${build_dir}/percolator/percolator*.rpm ${release_dir}
