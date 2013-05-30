#!/bin/bash

# execute the build scripts on a virtualbox host

pass="my_pass"
host="fedora"

ssh ${host} 'echo ' ${pass} ' | sudo -Sv && bash -s' < mingw64_build.sh
ssh ${host} 'echo ' ${pass} ' | sudo -Sv && bash -s' < rpm_build.sh

