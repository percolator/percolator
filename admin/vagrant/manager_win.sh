#!/bin/bash

#--------------------------------------------------------

usage()
{
    cat << EOF
This script uses vagrant to start a virtual machine and 
builds percolator with converters packages by executing 
a builder file.

usage: $0 
                  [[-h]] [[-a]]
                  [[-b branch]]|[[-s sourc_directory]]
                  [[-r release_directory]]
                  -p ubuntu|fedora|w32|w64|nw32|nw64

If no branch and source_directory is provided, the source
code from which the sourcecode is checked out from will be used.
Make sure that Vagrant and VirtualBox are up to date.
  -h     prints this help page
  -a     keeps the vagrant box alive (i.e. do not call vagrant destroy)
EOF
}
#--------------------------------------------------------
#input options management
#--------------------------------------------------------
#Default values:

current_path=$PWD;
script_dir=$(dirname ${BASH_SOURCE});
cd ${script_dir};
builder_adr="../builders/"
# boxes, might be overridden later on
vagbox_name="fedora18"
vagbox_url="http://www.nada.kth.se/~alinar/fedora18.box"

while getopts “hab:s:r:p:” OPTION; do
    case $OPTION in
        h)  usage;exit 1;;
        a)  alive="1";;
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
                nw32) post="nativew32"
                      batfile=true
                      vagbox_name="win7vs12"
                      vagbox_url="~/VagrantWin7/win7vs12.box"
                      ;;
                nw64) post="nativew64"
                      batfile=true
                      vagbox_name="win7vs12"
                      vagbox_url="~/VagrantWin7/win7vs12.box"
                      ;;
                *)
                     if [[ $OPTARG == *,* ]]; then
                         arr=$(echo $OPTARG | tr "," "\n");
                         multi_platform="1";
                     else
                         echo "Platform $OPTARG is undefined."
                         exit 1
                     fi;;
            esac;;
        \?)  echo "Invalid option: -${OPTARG}" >&2;;
    esac
done
if [[ -z $batfile ]]; then
  builder="${post}_build.sh";
else
  builder="${post}_build.bat";
fi

######
if [[ ! -z $multi_platform ]]; then
  [[ ! -z $branch ]] && arg=$arg"-b $branch "; 
  [[ ! -z $src ]] && arg=$arg"-s $src "; 
  [[ ! -z $release ]] && arg=$arg"-r $release "; 
  for x in $arr; do
    call_arg=$arg"-p $x";
    echo "Recusive call: $0 $call_arg"
    bash $0 $call_arg
  done
  exit 0;
fi
if [[ -z $post ]]; then
    usage
    echo "Please select one or more platforms with -p option."
    exit 1
fi
######

echo "------------------------------------------------------------------------";
echo "About to start up a build note using the builder script ${builder}";
echo "------------------------------------------------------------------------";


tmp_dir="$(mktemp -d --tmpdir tmp_${post}_XXXX)"
mkdir -p ${tmp_dir}/src/percolator
if [[ -z $src ]]; then
  if [[ -z $branch ]]; then
    echo "Copying source code from ${script_dir} to ${tmp_dir}/src/percolator/"
    cp -R ${script_dir}/../../* ${tmp_dir}/src/percolator/
  else
    echo "Cloning source code using the branch ${branch}"
    git clone --branch ${branch} https://github.com/percolator/percolator.git ${tmp_dir}/src/percolator
  fi
else
  echo "Copying source code from user specified path ${src}"
  cp -R ${src}/* ${tmp_dir}/src/percolator
fi
######
if [[ -z $release ]]; then
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

if [[ -z $batfile ]]; then
#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "2048", "--cpus", "4"]
  end
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = "-s /vagrant/src -r /vagrant/"
  end
end
EOF
#-----------------end of Vagrantfile content--------
#  config.vm.provision :shell, :inline => "su vagrant -c 'bash /vagrant/${builder} /vagrant/src /vagrant/build_${post}'"
else
cat <<EOF > Vagrantfile
Vagrant.configure("2") do |config|

  # Configure base box parameters
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.vm.guest = :windows
  config.winrm.username = "IEUser"
  config.winrm.password = "Passw0rd!"
  config.windows.halt_timeout = 30
  
  # Port forward WinRM and RDP
  config.vm.network :forwarded_port, guest: 3389, host: 3389
  config.vm.network :forwarded_port, guest: 5985, host: 5985, id: "winrm", auto_correct: true
  
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "2048", "--cpus", "4"]
  end
  
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = '-s "C:\vagrant\src" -b "C:\vagrant\build" -r "C:\vagrant"'
  end
end
EOF
fi
#---------------------------------------------------------------------------------------
vagrant up

#---------------------------------------------------------------------------------------
# release:

echo "Copying ready made packages from ${tmp_dir} to ${release}" 

mkdir -p ${release};
cp -v ${tmp_dir}/per*.{rpm,deb,exe,dmg} ${release};
cp -v ${tmp_dir}/elude*.{rpm,deb,exe,dmg} ${release};


#cp -v ${tmp_dir}/build_${post}/percolator/{per*.rpm,per*.deb,per*.exe,per*.dmg} ${release};
#cp -v ${tmp_dir}/build_${post}/converters/{per*.rpm,per*.deb,per*.exe,per*.dmg} ${release};
#---------------------------------------------------------------------------------------

if [[ -z ${alive} ]]; then
  vagrant destroy -f
fi

#---------------------------------------------------------------------------------------
