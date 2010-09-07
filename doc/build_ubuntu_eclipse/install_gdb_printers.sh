#!/bin/sh
# Mattia Tomasoni - Percolator Project
# Script that installs support for human readable print of STL containers while debugging with GDB
# Parameters: none

src=$(pwd)

sudo apt-get install subversion

cd $src
svn co svn://gcc.gnu.org/svn/gcc/trunk/libstdc++-v3/python


if [ "${HOME}/.gdbinit" ]; then
  rm ${HOME}/.gdbinit
fi

echo "python" >> ${HOME}/.gdbinit
echo "import sys" >> ${HOME}/.gdbinit
echo "sys.path.insert(0, '$src/python')" >> ${HOME}/.gdbinit
echo "from libstdcxx.v6.printers import register_libstdcxx_printers" >> ${HOME}/.gdbinit
echo "register_libstdcxx_printers (None)" >> ${HOME}/.gdbinit
echo "end" >> ${HOME}/.gdbinit
