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
                  [[-b branch]]|[[-s source_directory]]
                  [[-r release_directory]]
                  -p ubuntu|centos|fedora|win32|win64|osx

If no branch and source_directory is provided, the source
code from which the sourcecode is checked out from will be used.
Make sure that Vagrant and VirtualBox are up to date.
  -h     prints this help page
  -a     keeps the vagrant box alive (i.e. do not call vagrant destroy)
EOF
}

#--------------------------------------------------------
# input options management
#--------------------------------------------------------

# Project specific variables
tmp_src_dir=src/percolator
git_url=https://github.com/percolator/percolator.git
package_prefixes=()
package_prefixes+=(percolator)

# Default values:
current_path=$PWD;
script_dir=$(dirname ${BASH_SOURCE});
cd ${script_dir};
builder_adr="../builders/"
# boxes, might be overridden later on
vagbox_name="box-cutter/fedora23"
vagbox_url=""

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
                    vagbox_name="ubuntu/trusty64"
                    vagbox_url=""
                    package_ext="deb"
                    ;;
                fedora)	
                    post="fedora64"
                    package_ext="rpm"
                    ;;
                centos) 
                    post="centos64"
                    vagbox_name="bento/centos-7.2"
                    #vagbox_name="bento/centos-6.7"
                    vagbox_url=""
                    package_ext="rpm"
                    ;;
                win32) 
                    post="nativew32"
                    batfile=true
                    vagbox_name="win10vs15"
                    vagbox_url="~/VagrantWin7/win10vs15.box"
                    package_ext="exe"
                    ;;
                win64) 
                    post="nativew64"
                    batfile=true
                    vagbox_name="win10vs15"
                    vagbox_url="~/VagrantWin7/win10vs15.box"
                    package_ext="exe"
                    ;;
                osx)
                    post="osx64"
                    package_ext="pkg"
                    vagbox_name="osx-sierra-0.3.1"
                    vagbox_url="https://vagrant-osx.nyc3.digitaloceanspaces.com/osx-sierra-0.3.1.box"
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
mkdir -p ${tmp_dir}/${tmp_src_dir}
if [[ -z $src ]]; then
  if [[ -z $branch ]]; then
    echo "Copying source code from ${script_dir} to ${tmp_dir}/${tmp_src_dir}/"
    cp -R ${script_dir}/../../* ${tmp_dir}/${tmp_src_dir}/
  else
    echo "Cloning source code using the branch ${branch}"
    git clone --branch ${branch} ${git_url} ${tmp_dir}/${tmp_src_dir}
  fi
else
  echo "Copying source code from user specified path ${src}"
  cp -R ${src}/* ${tmp_dir}/${tmp_src_dir}
fi
######
if [[ -z $release ]]; then
  mkdir ${tmp_dir}/${post}_release
  release="${tmp_dir}/${post}_release"
fi

echo "Executing build procedure using:"
echo "tmp_dir=${tmp_dir}" 
echo "post=${post}" 

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

if [[ "$post" == "osx64" ]]; then

# run this on the ubuntu host to support nfs shares:
# sudo apt-get install nfs-common nfs-kernel-server

# copy the PackageMaker app
git clone https://github.com/erdnuesse/build-tools.git ${tmp_dir}/build-tools
# install the PackageMaker app on the host
sed -i '1i sudo hdiutil attach /vagrant/build-tools/xcode44auxtools6938114a.dmg; sudo cp -r /Volumes/Auxiliary\\ Tools/PackageMaker.app /Applications/' ${tmp_dir}/${builder}
# the nfs share is owned by the user on the client which blocks access by the vagrant user on the host machine
# currently, there does not seem to be a way to map the user ids for the nfs share...
# instead, force hdiutil to use sudo to force access to the nfs share
sed -i '1i sudo printf "#!/bin/bash\\nsudo /usr/bin/hdiutil \\\"\\$@\\\"" > /usr/local/bin/hdiutil; sudo chmod +x /usr/local/bin/hdiutil' ${tmp_dir}/${builder}

#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.ssh.insert_key = false
  #config.ssh.password = "vagrant" # private key authentication does not work on this box
  config.vm.boot_timeout = 600
  config.vm.network "private_network", ip: "192.168.56.10"
  
  # Synced folder are not supported under Mac OS X
  # config.vm.synced_folder ".", "/vagrant", :disabled => true
  # Use NFS for the shared folder
  config.vm.synced_folder ".", "/vagrant",
    id: "vagrant-root",
    :nfs => true,
    :mount_options => ['nolock,vers=3,udp,noatime,actimeo=1,resvport'],
    :export_options => ['async,insecure,no_subtree_check,no_acl,no_root_squash']

  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "4096", "--cpus", "4"]
    # vb.gui = true # turn on for trouble shooting, e.g. if boot times out repeatedly
    # Fix "hfs mounted macintosh hd on device root_device" issue
    vb.customize ["modifyvm", :id, "--cpuidset", "1","000206a7","02100800","1fbae3bf","bfebfbff"]

    # Some more hacks for device recognition
    vb.customize ["setextradata", :id, "VBoxInternal/Devices/efi/0/Config/DmiSystemProduct", "MacBookPro11,3"]
    vb.customize ["setextradata", :id, "VBoxInternal/Devices/efi/0/Config/DmiSystemVersion", "1.0"]
    vb.customize ["setextradata", :id, "VBoxInternal/Devices/efi/0/Config/DmiBoardProduct", "Iloveapple"]
    vb.customize ["setextradata", :id, "VBoxInternal/Devices/smc/0/Config/DeviceKey", "ourhardworkbythesewordsguardedpleasedontsteal(c)AppleComputerInc"]
    #vb.customize ["setextradata", :id, "VBoxInternal/Devices/smc/0/Config/GetKeyFromRealSMC", "1"]
  end
  
  config.vm.provision :shell do |shell|
    shell.privileged = false
    shell.path = "${tmp_dir}/${builder}"
    shell.args = "-s /vagrant/src -b /vagrant/build -r /vagrant/"
  end
end
EOF

elif [[ -z $batfile ]]; then

vagbox_url_line=""
if [[ -n ${vagbox_url} ]]; then
  vagbox_url_line="config.vm.box_url = \"${vagbox_url}\""
fi

#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "${vagbox_name}"
  ${vagbox_url_line}
  config.ssh.insert_key = false
  config.vm.boot_timeout = 600
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "4096", "--cpus", "4"]
    # vb.gui = true # turn on for trouble shooting, e.g. if boot times out repeatedly
  end
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = "-s /vagrant/src -b /vagrant/build -r /vagrant/"
  end
end
EOF
#-----------------end of Vagrantfile content--------
#  config.vm.provision :shell, :inline => "su vagrant -c 'bash /vagrant/${builder} /vagrant/src /vagrant/build_${post}'"
else
#-----------------mount_vagrant_drive.bat content---------------

# Qt build currently fails on a shared drive, but works on a mapped network drive
sed -i '1i net use z: \\\\vboxsrv\\vagrant' ${tmp_dir}/${builder}  

cat <<EOF > Vagrantfile
Vagrant.configure("2") do |config|

  # Configure base box parameters
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.vm.guest = :windows
  config.windows.halt_timeout = 30
  config.vm.boot_timeout = 1200

  # Port forward WinRM and RDP
  config.vm.communicator = "winrm"
  # config.vm.network "forwarded_port" , host: 33390 , guest: 3389 # allows remote desktop with "vagrant rdp"
  
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "8192", "--cpus", "4"]
    # vb.gui = true # turn on for trouble shooting, e.g. if boot times out repeatedly
  end
  
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    #shell.args = '-s "C:\vagrant\src" -b "C:\vagrant\build" -r "C:\vagrant"'
    shell.args = '-s "C:\vagrant\src" -b "Z:\build" -r "C:\vagrant"'
  end
end
EOF
fi
#---------------------------------------------------------------------------------------
vagrant up

if [[ $? -eq 0 ]]; then
  echo "Building of binaries succeeded"
else
  echo "Building of binaries failed"
  alive="1"
fi

#---------------------------------------------------------------------------------------
# release:

echo "Copying ready made packages from ${tmp_dir} to ${release}" 

mkdir -p ${release};
for package_prefix in ${package_prefixes[@]}; do
  cp -v ${tmp_dir}/${package_prefix}*.${package_ext} ${release};
done

#---------------------------------------------------------------------------------------

if [[ -z ${alive} ]]; then
  vagrant destroy -f
else
  echo "-a option set or encountered error: keeping the VM alive, remember to close and delete the VM manually."
  exit 1
fi

#---------------------------------------------------------------------------------------
