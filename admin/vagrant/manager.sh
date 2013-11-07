#!/bin/bash

#--------------------------------------------------------

usage()
{
    cat << EOF
This script uses vagrant to start a virtual machine and 
builds percolator with converters packages by executing 
a builder file.

usage: $0 
                  [-h] (print this message)
                  [-b branch]|[-s sourc_directory]
                  [-r release_directory]
                  -p ubuntu|fedora|w32|w64

If no branch and source_directory is provided, the same source
code will be used.
Make sure that Vagrant and VirtualBox are up to date.
EOF
}
#--------------------------------------------------------
#input options management
#--------------------------------------------------------
#Default values:

current_path=$PWD;
cd $(dirname ${BASH_SOURCE});
source_path=$PWD;
multi_platform=0
builder_adr="../builders/"
# boxes, might be overridden later on
vagbox_name="fedora18"
vagbox_url="http://www.nada.kth.se/~alinar/fedora18.box"

while getopts “hb:s:r:p:” OPTION; do
    case $OPTION in
        h)  usage;exit 1;;
        b)  branch=$OPTARG;;
        s)  src=$OPTARG;;
        r)  release=$OPTARG;;
        p)  case $OPTARG in
               	ubuntu)
                    post="ubuntu64"
                    vagbox_name="precise64"
                    vagbox_url="http://files.vagrantup.com/precise64.box"
                    ;;
                fedora)	post="fedora64";;
                w64) post="mingw64";;
                w32) post="mingw32";;
                *)
                     if [[ $OPTARG == *,* ]]; then
                         arr=$(echo $OPTARG | tr "," "\n")
                         multi_platform=1
                     else
                         echo "Platform $OPTARG is undefined."
                         exit 1
                     fi;;
            esac;;
        \?)  echo "Invalid option: -${OPTARG}" >&2;;
    esac
done
builder="${post}_build.sh"


######
if [ $multi_platform == 1 ]; then
  for x in $arr; do
    arg=""
    if [ ! -z $branch ]; then
      arg=$arg"-b $branch " 
    fi
    if [ ! -z $src ]; then
      arg=$arg"-s $src " 
    fi
    if [ ! -z $release ]; then 
      arg=$arg"-r $release " 
    fi
    arg=$arg"-p $x"
    # recusive call to this script, but specifying one platform.
    bash $0 $arg
  done
  exit 0;
fi
######
if [ -z $post ]
then
    usage
    echo "Please select one or more platforms with -p option."
    exit 1
fi
######
tmp_dir="$(mktemp -d --tmpdir tmp_${post}_XXXX)"
mkdir ${tmp_dir}/src ${tmp_dir}/src/percolator
if [ -z $src ]; then
  if [ -z $branch ]; then
    echo "Copying source code from ${source_path}"
    cp -R ${source_path}/../../* ${tmp_dir}/src/percolator
  else
    echo "Cloning source code using the branch ${branch}"
    git clone --branch ${branch} https://github.com/percolator/percolator.git ${tmp_dir}/src/percolator
  fi
else
  echo "Copying source code from user specified path ${src}"
  cp -R ${src}/* ${tmp_dir}/src/percolator
fi
######
if [ -z $release ]; then
  mkdir ${tmp_dir}/${post}_release
  release="${tmp_dir}/${post}_release"
fi

echo "Executing build procedure using:"
echo "tmp_dir=${tmp_dir}" 
echo "post=${post}" 
echo "tmp_dir=${tmp_dir}" 


#--------------------------------------------------------
#########################################################
#--------------------------------------------------------
#copy builder:
cp ${builder_adr}${builder} ${tmp_dir};
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# making the Vagrantfile:
cd ${tmp_dir};
touch Vagrantfile;
#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = "-s /vagrant/src -r /vagrant/"
  end
end
EOF
#-----------------end of Vagrantfile content--------
#  config.vm.provision :shell, :inline => "su vagrant -c 'bash /vagrant/${builder} /vagrant/src /vagrant/build_${post}'"

#---------------------------------------------------------------------------------------
vagrant up

#---------------------------------------------------------------------------------------
# release:

echo "Copying ready made packages from ${tmp_dir} to ${release}" 

mkdir -p ${release};
cp -v ${tmp_dir}/per*.{rpm,deb,exe,dmg} ${release};

#cp -v ${tmp_dir}/build_${post}/percolator/{per*.rpm,per*.deb,per*.exe,per*.dmg} ${release};
#cp -v ${tmp_dir}/build_${post}/converters/{per*.rpm,per*.deb,per*.exe,per*.dmg} ${release};
#---------------------------------------------------------------------------------------

vagrant destroy -f
#---------------------------------------------------------------------------------------
