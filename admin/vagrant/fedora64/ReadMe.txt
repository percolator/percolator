Prerequisite:
. VirtualBox
. git


Instructions:
1. First you need to import the fedora64 virtual machine (VM) to VirtualBOX. To do so, just run '. import_new_vm.sh'. it automatically downloads 'fedora64.ova' (VM file) and imports it in VirtualBox. Remember to remove any VM in Virtualbox with the name 'fedora64' before running this file. This step is needed only for the first time.

2. in build.sh change the release value to desierd subfolder. it is recommended to leave the other values as they are.

3. run '. build.sh'


