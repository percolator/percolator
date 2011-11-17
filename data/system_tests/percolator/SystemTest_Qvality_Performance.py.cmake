# Mattia Tomasoni - Percolator Project
# Script that tests the performances of qvality
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
success = True

print "QVALITY PERFORMANCE"

# the output line containing "Selecting pi_0" is extracted and if its value is 
# outside of (0.86, 0.90) an error is reported
print "(*): checking selected pi_0..."
processFile = os.popen("grep \"Selecting pi_0\" " + 
  "/tmp/qvalityOutput.txt")
output = processFile.read()
extracted = float(output[15:20])
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0=" + str(extracted) + " is outside of desired range (0.86, 0.90)"
  print "check /tmp/qvalityOutput.txt for details" 
  success = False

# the number of lines of stdout (after the line beginning with "Score") until 
# q-value < 0.01 are counted and an error is reported if their number is greater
# than 755+/-5%
print "(*): checking values..."
processFile = open("/tmp/qvalityOutput.txt")
line = processFile.readline()
finished = False
while (not finished): # reading line by line, looking for "Score"
  if line[0:5] == "Score":
    line = processFile.readline()
    finished = True
  else:
    line = processFile.readline()
countLines = 0
finished = False
while (not finished): # counting lines
  extracted = ""
  tabsNumber = 0
  for i in range(0, len(line)):
    if line[i] == "\t":
      tabsNumber = tabsNumber + 1
    if tabsNumber == 2:
      extracted += line[i] # right column extracted
  extracted = float(extracted)
  if extracted > 0.01:
    finished = True
  else:
    countLines = countLines + 1
    line = processFile.readline()
if countLines < 717 or countLines > 793:
  print "...TEST FAILED: number of peptides=" + str(countLines) + " outside of desired range (717, 793)"
  print "check /tmp/qvalityOutput.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)
