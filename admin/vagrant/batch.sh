#!/bin/bash
export LANG=en_US.UTF-8
export LC_ALL="en_US.UTF-8"
trap 'echo "Batch script killed"; exit 1' INT TERM

release_dir=${HOME}/release

platforms=()
platforms+=(osx)
platforms+=(ubuntu)
platforms+=(centos)
platforms+=(fedora)
platforms+=(win64)
platforms+=(win32)

for platform in ${platforms[@]}; do
  echo "Building $platform binaries"
  echo "  Check $(pwd)/${platform}_output.txt for progress"
  ./manager.sh -p $platform -r ${release_dir}/$platform > ${platform}_output.txt 2>&1
  if [[ $? -eq 0 ]]; then
    echo "Building of ${platform} binaries succeeded"
  else
    echo "Building of ${platform} binaries failed"
  fi
done

# for the osx build, root priviliges are necessary but we can skip the password prompt:
# Vagrant NFS access https://www.vagrantup.com/docs/synced-folders/nfs.html#root-privilege-requirement
# add the following lines with $sudo visudo
#Cmnd_Alias VAGRANT_EXPORTS_CHOWN = /bin/chown 0\:0 /tmp/*
#Cmnd_Alias VAGRANT_EXPORTS_MV = /bin/mv -f /tmp/* /etc/exports
#Cmnd_Alias VAGRANT_NFSD_CHECK = /etc/init.d/nfs-kernel-server status
#Cmnd_Alias VAGRANT_NFSD_START = /bin/systemctl start nfs-server.service
#Cmnd_Alias VAGRANT_NFSD_APPLY = /usr/sbin/exportfs -ar
#%sudo ALL=(root) NOPASSWD: VAGRANT_EXPORTS_CHOWN, VAGRANT_EXPORTS_MV, VAGRANT_NFSD_CHECK, VAGRANT_NFSD_START, VAGRANT_NFSD_APPLY
