## Percolator
[![Build status](https://github.com/percolator/percolator/actions/workflows/build.yml/badge.svg)](https://github.com/percolator/percolator/actions/workflows/build.yml) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md)

This repository holds the source code for percolator and format
converters software packages; a software for postprocessing of shotgun
proteomics data.

The source code is hosted at github, binary downloads can be found in the release repository:
https://github.com/percolator/percolator/releases

The building procedure is cmake-based. Helper scripts for building on
different platforms can be found in: admin/builders/

    cmake -S . -B build
    cmake --build build -j 4
    cd build && ctest

Requirements: CMake 3.15 or later and a C++17 compiler. The oldest
supported GCC is 9; from GCC 11 on, reading and writing numbers use the
faster std::from_chars and std::to_chars. The continuous integration builds
with GCC 11 and 13 on Linux (x86-64 and arm64), AppleClang on macOS 15 (the
binaries run on macOS 12 and later) and Visual Studio 2022 on Windows; its
job with GCC 9 does not block yet.
Results can differ in a printed digit between platforms, since compilers and
math libraries round differently; the Build workflow lists the differences
(job platform-differences).
