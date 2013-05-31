#!/bin/bash

# execute the build scripts on a virtualbox host

pass="my_pass" # the pass word for the sudo rights on the VBox
host="fedora" # ip-address of the VBox

ssh ${host} 'echo ' ${pass} ' | sudo -Sv && bash -s' < rpm_build.sh
ssh ${host} 'echo ' ${pass} ' | sudo -Sv && bash -s' < mingw64_build.sh

