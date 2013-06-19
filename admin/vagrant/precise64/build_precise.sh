#!/bin/bash

branch="branch-2-05";
release="$HOME/release";
#---------------------------------------------------------------------------------------
tmp_dir="/tmp/tmp";
builder="precise64_build.sh";
precise_box_url="http://files.vagrantup.com/precise64.box";
#---------------------------------------------------------------------------------------
# making directories and copy builder:
rm -rf ${tmp_dir}; mkdir ${tmp_dir};
mkdir ${tmp_dir}/src;
mkdir -p ${release};
cp ${builder} ${tmp_dir};
#---------------------------------------------------------------------------------------
# clone :
git clone --branch ${branch} https://github.com/percolator/percolator.git ${tmp_dir}/src/percolator;

#---------------------------------------------------------------------------------------
# making the Vagrantfile:
cd ${tmp_dir};
touch Vagrantfile;
#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "precise64"
  config.vm.box_url = "$precise_box_url"
  config.vm.provision :shell, :path => "${builder}"
end
EOF
#-----------------end of Vagrantfile content--------

#---------------------------------------------------------------------------------------
vagrant up

#---------------------------------------------------------------------------------------
# release:
cp ${tmp_dir}/build/percolator/per*.deb ${release};
cp ${tmp_dir}/build/converters/per*.deb ${release};

#---------------------------------------------------------------------------------------
vagrant destroy -f
#---------------------------------------------------------------------------------------
cd ${release};

