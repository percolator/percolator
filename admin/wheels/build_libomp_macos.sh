#!/usr/bin/env bash
# Builds a static LLVM OpenMP runtime (libomp.a) for the macOS deployment
# target of the Percolator binaries, for the wheel build (pyproject.toml,
# [tool.cibuildwheel.macos]).
#
# Homebrew's libomp is built for the macOS version of the build machine.
# Linked into binaries with an older deployment target (12.0), it makes the
# linker warn about the newer object files, and nothing checks that it runs
# on the older system. This script compiles the runtime from the LLVM
# sources for MACOSX_DEPLOYMENT_TARGET (default 12.0) and the architecture
# of the machine. The two source archives are checked against the SHA-256
# sums below before anything is unpacked.
#
# usage: admin/wheels/build_libomp_macos.sh <install-prefix>
#   Gives <install-prefix>/lib/libomp.a and <install-prefix>/include/omp.h.
#   Pass them to the Percolator build with
#     -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I<install-prefix>/include"
#     -DOpenMP_CXX_LIB_NAMES=omp
#     -DOpenMP_omp_LIBRARY=<install-prefix>/lib/libomp.a
#   Needs curl, cmake and a C++ compiler. The script holds nothing specific
#   to macOS except the default deployment target; with clang on Linux it
#   builds the same library, which is how it was tried.
set -euo pipefail
if [ $# -ne 1 ]; then
  sed -n '2,22p' "$0" >&2
  exit 2
fi
mkdir -p "$1"
PREFIX=$(cd "$1" && pwd)

LLVM_VERSION=19.1.7
OPENMP_SHA256=bd7e6901ab086fd268750363017935fd4a717c153dad3c2aab86cb0140d9e3fe
CMAKE_SHA256=11c5a28f90053b0c43d0dec3d0ad579347fc277199c005206b963c19aae514e3
URL=https://github.com/llvm/llvm-project/releases/download/llvmorg-$LLVM_VERSION

export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-12.0}

sha256() {
  if command -v sha256sum > /dev/null; then
    sha256sum "$1" | cut -d ' ' -f 1
  else
    shasum -a 256 "$1" | cut -d ' ' -f 1
  fi
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# the runtime's sources and the CMake modules they share with LLVM, which
# the build expects in ../cmake
fetch() {  # component, expected sum
  local archive=$1-$LLVM_VERSION.src.tar.xz
  curl -fsSL --retry 3 -o "$archive" "$URL/$archive"
  local sum
  sum=$(sha256 "$archive")
  if [ "$sum" != "$2" ]; then
    echo "error: $archive has SHA-256 $sum, expected $2" >&2
    exit 1
  fi
  mkdir "$1"
  tar -xJf "$archive" -C "$1" --strip-components 1
}
fetch openmp "$OPENMP_SHA256"
fetch cmake "$CMAKE_SHA256"

cmake -S openmp -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DLIBOMP_ENABLE_SHARED=OFF \
  -DLIBOMP_INSTALL_ALIASES=OFF \
  -DLIBOMP_OMPT_SUPPORT=OFF \
  -DOPENMP_ENABLE_LIBOMPTARGET=OFF \
  -DOPENMP_ENABLE_OMPT_TOOLS=OFF \
  -DOPENMP_ENABLE_TESTING=OFF
cmake --build build -j 4
cmake --install build

ls -l "$PREFIX/lib/libomp.a" "$PREFIX/include/omp.h"
