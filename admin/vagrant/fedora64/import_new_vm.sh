#!/bin/bash

#the URL address of the VM
addr="http://www.nada.kth.se/~alinar/fedora64.ova";

#check for existance of fedora64 on virtualbox
if (vboxmanage list vms | grep fedora64) &> /dev/null;
then echo "fedora64 already exists. Please remove it before importing a new vm.";
return 1;
fi;
#
echo "Please make sure that there is no VM in VirtualBox with name 'fedora64'."

#download the VM
wget -P /tmp $addr;

#import the VM
VBoxManage import /tmp/fedora64.ova;

#remove the downloaded file
rm /tmp/fedora64.ova
