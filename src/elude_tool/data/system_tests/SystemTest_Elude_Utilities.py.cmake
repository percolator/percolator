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
  print(message)
  print("...TEST FAILED")
  
# check the existence of a list of files 
# the test_name is used for printing 
def checkFilesExistence(test_name, files):
  for f in files:
    if not os.path.isfile(f) :
      failTest("{0}, unable to locate file {1}".format(test_name, f))
      return False
  return True  
      
# clean up by deleting a list of files 
def cleanUp(files): 
  for f in files: 
    try:
      os.remove(f) 
    except:
      pass

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
    return (None, None, None)
     
# check that the index file includes the given symbols at the given lines 
def checkIndex(test_name, index_file, indices, values):
  lines = loadFile(index_file)
  idx = map(lambda line: line.split(" : ")[0].strip(), lines)
  for index, val in zip(indices, values):
    if (idx[index] != val):    
      failTest(test_name + ", incorrect symbols in the index")
      return False 
  return True
    
# check content of the output file 
def checkOutputFile(test_name, out_file, no_lines = -1, indices = None, 
                    epeps = None, eprts = None, eorts = None):
  lines = loadFile(out_file)
  sorted(lines)
  if (no_lines != -1 and len(lines) != no_lines): 
    failTest(test_name + ", incorrect number of peptides in the output file")
    return False
  
  if (indices != None):
    peps = map(lambda i: lines[i].split("\t")[0], indices)
    prts = map(lambda i: float(lines[i].split("\t")[1]), indices)
    if (eorts != None): 
      orts = map(lambda i: float(lines[i].split("\t")[2]), indices) 
    
    for i in range(len(indices)):
      if peps[i] != epeps[i]:     
        failTest(test_name + ", incorrect peptides in the output file")
        return False
      if (abs(prts[i] - eprts[i] > eprts[i]*0.15)):       
        failTest(test_name + ", incorrect predictions in the output file")
        return False
      if (eorts != None) and (abs(orts[i] - eorts[i] > 0.1)):               
        failTest(test_name + ", incorrect observed rts in the output file")
        return False
  return True
