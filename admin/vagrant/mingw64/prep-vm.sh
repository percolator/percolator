#!/bin/bash
vm_name="fedora64";
share_dir=$1;

#poweroff fedora64 vm.
echo "Shutting down fedora64 VM ...";
(VBoxManage controlvm fedora64 poweroff) &> /dev/null;
 
#checking no VM is running
#if [ $( VBoxManage list runningvms | wc -l ) -ne 0 ];
#then echo "Please shutdown ALL the VMs in VirtualBox and try again.";return 1;fi;

#set port forwarding for ssh
echo "set port forwarding 22 -> 2222 ...";
(VBoxManage modifyvm $vm_name --natpf1 "guestssh,tcp,,2222,,22") &> /dev/null;

#adding the share directories
echo "Adding shared folders to VM ...";
(VBoxManage sharedfolder add $vm_name --name "share-dir" --hostpath $share_dir) &> /dev/null;

#start running fedora64 VM
(VBoxManage startvm fedora64 --type headless) &> /dev/null;

#using ssh with the key to mount the shared forlder into /shared in the VM
echo "Mounting the shared folders in VM. This may take a few minutes......";
(ssh -i insecure_private_key -l root -p 2222 localhost 'bash -s' < vm-mount.sh) &> /dev/null;




