description:
this script makes cross-compilation using mingw64 in a fedora64 vm and releases WIN32 exe packages of percolator and converters.

Prerequisite:
. VirtualBox
. git


Instructions:
1. First you need to import the fedora64 virtual machine (VM) to VirtualBOX. To do so, just run '. import_new_vm.sh'. it automatically downloads 'fedora64.ova' (VM file) and imports it in VirtualBox. Remember to remove any VM in Virtualbox with the name 'fedora64' before running this file. This step is needed only for the first time.

2. In build_mingw.sh change the release value to desierd subfolder and change branch to the desierd branch. it is recommended to leave the other values as they are.

3. run '. build_mingw.sh'


