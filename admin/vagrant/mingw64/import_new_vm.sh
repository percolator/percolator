#!/bin/bash

#the URL address of the VM
addr="http://www.nada.kth.se/~alinar/fedora64.ova";

#
echo "Please make sure that there is no VM in VirtualBox with name 'fedora64'."

#download the VM
wget -P /tmp $addr;

#import the VM
VBoxManage import /tmp/fedora64.ova;

#remove the downloaded file
rm /tmp/fedora64.ova
