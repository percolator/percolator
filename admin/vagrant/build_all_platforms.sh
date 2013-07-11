#!/bin/bash
branch="branch-2-05";
release="$HOME/release";
#-------------------------------------------------------------
tmp_dir="$(mktemp -d --tmpdir master_tmp_XXXX)";
mkdir ${tmp_dir}/src;
mkdir -p ${release};
#-------------------------------------------------------------
Ppath="$(dirname ${BASH_SOURCE})";
cd ${Ppath};
Spath=$PWD;
#-------------------------------------------------------------
# clone :
git clone --branch ${branch} https://github.com/percolator/percolator.git ${tmp_dir}/src/percolator;
#-------------------------------------------------------------
cd ${Spath}/precise64;bash build_precise.sh ${tmp_dir} ${release}/precise;
cd ${Spath}/fedora;bash build_fedora.sh ${tmp_dir} ${release}/fedora;
cd ${Spath}/mingw;bash build_win.sh ${tmp_dir} ${release}/mingw;
#-------------------------------------------------------------

