#!/bin/bash

release=$HOME;
branch="branch-2-05";
tmpdir="/tmp/tmp";

#insecure_private_key download and permission adjustments
if [ ! -f insecure_private_key ];
then wget http://www.nada.kth.se/~alinar/insecure_private_key;
fi;
chmod 600 insecure_private_key;

#clone into source----------------------------------------
echo "Cloning Percolator branch:$branch.....";
rm -rf $tmpdir;mkdir $tmpdir;
git clone https://github.com/percolator/percolator.git $tmpdir/percolator --branch $branch;

#start and config the vm-----------------------------------
. prep-vm.sh $tmpdir;

#running rpm_build in vm-----------------------------------
ssh -i insecure_private_key -l root -p 2222 localhost 'bash -s' < mingw64_build.sh

#release--------------------------------------------------
mkdir $release;
cp $tmpdir/build/percolator/per*.exe $release;
cp $tmpdir/build/converters/per*.exe $release;
#poweroff fedora64 vm--------------------------------------
echo "Shutting down fedora64 VM ...";
(VBoxManage controlvm fedora64 poweroff) &> /dev/null;
