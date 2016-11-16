#!/bin/bash
export LANG=en_US.UTF-8
export LC_ALL="en_US.UTF-8"

echo "Building Ubuntu binaries"
./manager.sh -p ubuntu > ubuntu_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of Ubuntu binaries succeeded"
fi

echo "Building Fedora binaries"
./manager.sh -p fedora > fedora_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of Fedora binaries succeeded"
fi

echo "Building CentOS binaries"
./manager.sh -p centos > centos_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of CentOS binaries succeeded"
fi

echo "Building native 64-bit Windows binaries"
./manager.sh -p nw64 > nw64_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of native 64-bin Windows binaries succeeded"
fi

echo "Building native 32-bit Windows binaries"
./manager.sh -p nw32 > nw32_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of 32-bit native Windows binaries succeeded"
fi

#./manager.sh -p w32 > w32_output.txt 2>&1 # MINGW32
#./manager.sh -p w64 > w64_output.txt 2>&1 # MINGW64
