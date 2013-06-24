#!/bin/bash

branch="branch-2-05";
release="$HOME/release";


#-------------------------------------------------------------
tmp_dir="$(mktemp -d --tmpdir master_tmp_XXXX)";
mkdir ${tmp_dir}/src;
#-------------------------------------------------------------
current_path=$PWD;
#-------------------------------------------------------------
# clone :
git clone --branch ${branch} https://github.com/percolator/percolator.git ${tmp_dir}/src/percolator;
#-------------------------------------------------------------
cd $current_path/precise64;bash build_precise.sh ${tmp_dir} ${release};
cd $current_path/fedora;bash build_fedora.sh ${tmp_dir} ${release};
cd $current_path/mingw;bash build_win.sh ${tmp_dir} ${release};
#-------------------------------------------------------------

