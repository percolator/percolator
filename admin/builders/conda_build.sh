#!/bin/bash
# Builds percolator and qvality, dynamically linked, for a conda package
# (bioconda) on Linux and macOS. The recipe's build.sh calls it:
#
#   bash admin/builders/conda_build.sh
#
# It uses PREFIX, CPU_COUNT and CMAKE_ARGS from conda-build. CMAKE_ARGS
# carries the settings of the conda compilers (sysroot, macOS deployment
# target, cross-compilation). To try it outside conda-build, set PREFIX to
# an installation directory.
set -euxo pipefail

src_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
build_dir="${src_dir}/build-conda"
: "${PREFIX:?PREFIX must name the installation directory}"

# Start from an empty build directory. conda-build copies the source tree,
# including a build-conda left by an earlier run of this script, and CMake
# refuses a cache made for another source directory.
rm -rf "${build_dir}"

# CMP0125=NEW: on Linux, CMAKE_ARGS names the conda binutils without a path
# (-DCMAKE_AR=x86_64-conda-linux-gnu-ar). With the old behaviour, which
# cmake_minimum_required(3.15) leaves in place, CMake makes such a name a
# path in the current directory, and the link of libperclibrary.a fails
# with "Error running link command: no such file or directory".
# STATIC=OFF: link OpenMP from the environment dynamically.
# HOMEBREW_LIBOMP=OFF: never link to a libomp outside the environment.
# NO_TESTS=ON: no tests in the package build; the unit tests would download
# GoogleTest.
# shellcheck disable=SC2086  # CMAKE_ARGS holds several arguments
cmake -S "${src_dir}" -B "${build_dir}" ${CMAKE_ARGS:-} \
  -DCMAKE_POLICY_DEFAULT_CMP0125=NEW \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DSTATIC=OFF \
  -DHOMEBREW_LIBOMP=OFF \
  -DNO_TESTS=ON
cmake --build "${build_dir}" --parallel "${CPU_COUNT:-1}"
cmake --install "${build_dir}"
