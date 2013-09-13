#!/bin/bash


sudo groupadd admin

pass=$(perl -e 'print crypt($ARGV[0], "password")' vagrant)
sudo useradd -m -p $pass -G admin vagrant 

sudo sed -i 's/Defaults.*requiretty/ /g' /etc/sudoers
sudo bash -c "echo '%admin ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers"

sudo yum install ssh

sudo mkdir -p /home/vagrant/.ssh
sudo bash -c "echo ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA6NF8iallvQVp22WDkTkyrtvp9eWW6A8YVr+kz4TjGYe7gHzIw+niNltGEFHzD8+v1I2YJ6oXevct1YeS0o9HZyN1Q9qgCgzUFtdOKLv6IedplqoPkcmF0aYet2PkEDo3MlTBckFXPITAMzF8dJSIFo9D8HfdOV0IAdx4O7PtixWKn5y2hMNG0zQPyUecp4pzC6kivAIhyfHilFR61RGL+GPXQ2MWZWFYbAGjyiYJnAmCP3NOTd0jMZEnDkbUvxhMmBYSdETk1rRgm+R4LOzFUGaHqHDLKLX+FIPKcF96hrucXzcWyLbIbEgE98OHlnVYCzRdK8jlqm8tehUc9c9WhQ== vagrant insecure public key >> /home/vagrant/.ssh/authorized_keys"

sudo chmod 600 /home/vagrant/.ssh/authorized_keys 
