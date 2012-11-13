#!/usr/bin/python
# @Created by L. Moruz
# November 2nd, 2010
# this script includes small function used in the tests
import os 
import sys
import subprocess

pathToBinaries = os.path.join("@pathToBinaries@", "gtest_unit_elude")

def main():

  ret = subprocess.call(pathToBinaries)
  if (ret != 0):
    exit(1)

if __name__ == '__main__':
  main()
  
