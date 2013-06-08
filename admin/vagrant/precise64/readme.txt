prerequisites:

.vagrant
.virtualbox
.ssh (for windows a software that contains ssh e.g. git can be installed and its bin folder should be added to path folders of windows)

-----------------------------------------------------------------------------
Instuctions:

.create a folder X
.copy 'vagrantfile' and 'script.sh' to the folder X.
.change branch value in script.sh to desired branch you want to be built.
.$ vagrant up
.a folder named 'build' is created containing built files. downloaded source folder will be there too.
.use '$vagrant halt' to turn off the VM and preserving installed packages OR use '$vagrant destraoy' to completely reset the VM.

-----------------------------------------------------------------------------

