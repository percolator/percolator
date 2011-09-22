#!/usr/bin/python
# @Created by tomasoni
# Feb 15th, 2011
# this script includes small function used in the tests
import os 
import sys
import subprocess


def main():
  pathToBinaries = "@pathToBinaries@"
  currPath = os.getcwd()
  os.chdir(pathToBinaries)
  ret = subprocess.call("./gtest_unit")
  os.chdir(currPath)
  if (ret != 0):
    exit(1)

if __name__ == '__main__':
  main()
  
