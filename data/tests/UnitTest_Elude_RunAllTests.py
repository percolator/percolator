#!/usr/bin/python
# @Created by L. Moruz
# November 2nd, 2010
# this script includes small function used in the tests
import os 
import sys
import subprocess

def main():
  path = os.path.join(os.path.dirname(sys.argv[0]), "../../")
  curr_path = os.getcwd()
  os.chdir(path)
  ret = subprocess.call("gtest_unit")
  os.chdir(curr_path)
  if (ret != 0):
    exit(1)

if __name__ == '__main__':
  main()
  