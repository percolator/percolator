#!/usr/bin/python
# @Created by L. Moruz
# October 15th, 2010 
# this script includes small function used in the tests
import os 

########################## Useful functions ##########################
# deletes the contents of a directory 
def removeDir(dirname):
  files = os.listdir(dirname)
  for f in files:
    os.remove(os.path.join(dirname, f))
  os.rmdir(dirname) 

# return the content of a file as a list of lines 
def loadFile(filename):
  handle = open(filename, "r")
  lines = handle.readlines()
  handle.close() 
  return lines
  
# check if the file filename includes no_lines lines and that the content of lines line_numbers is line_contents
def testFileContent(filename, no_lines, line_numbers, line_contents):
  lines = loadFile(filename)
  
  # check the number of lines 
  if (len(lines) != no_lines):
    return False
  
  # check the content of the lines 
  for l, lc in zip(line_numbers, line_contents):
    if (lines[l] != lc):
      return False
    
  return True

# fail a test with a message 
def failTest(message):
  print message 
  print "...TEST FAILED"
  exit(1)
  
# check the existence of a list of files 
# the test_name is used for printing 
def checkFilesExistence(test_name, files):
  for f in files:
    if not os.path.isfile(f) :
      failTest("{0}, unable to locate file {1}".format(test_name, f))
      exit(1) 
      
# clea up by deleting a list of files 
def cleanUp(files): 
  map(os.remove, files) 

# check the performance measures from a log file
# values includes the expected pcorr, scorr and delta
def checkPerformance(test_name, filename, values = None): 
  lines = loadFile(filename)
  pcorr = 0.0
  scorr = 0.0
  delta = 1000
  for l in lines: 
    if l.strip().startswith("Pearson's"): 
      pcorr = float(l.split("=")[1])
    if l.strip().startswith("Spearman's"): 
      scorr = float(l.split("=")[1])
    if l.strip().startswith("Delta_t"): 
      delta = float(l.split("=")[1])
 
  if (values == None):
    return (pcorr, scorr, delta)
  (p, s, d) = values 
  if pcorr < p - 0.15 or scorr < s - 0.15 or delta > d + (d * 0.15): 
    failTest(test_name + ", incorrect performance figures")
     