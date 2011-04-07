# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
success = True

print "PERCOLATOR PERFORMANCE"

# the output line containing "New pi_0 estimate" is extracted and if its value is 
# outside of 622+/-5% an error is reported
print "(*): checking new pi_0 estimate in Percolator's output..."
processFile = os.popen("grep \"New pi_0 estimate\" " + 
  "/tmp/PERCOLATOR_no_options.txt")
output = processFile.read()
extracted = float(output[39:42])
if extracted < 590.9 or extracted > 653.1: 
  print "...TEST FAILED: new pi_0 estimate=" + str(extracted) + " is outside of desired range (590.9, 653.1)"
  print "check /tmp/PERCOLATOR_no_options.txt for details" 
  success = False

# the output line containing "Selecting pi_0" is extracted and if its value is 
# outside of (0.86, 0.90) an error is reported
print "(*): checking selected new pi_0 estimate in Percolator's output..."
processFile = os.popen("grep \"Selecting pi_0\" " +
  "/tmp/PERCOLATOR_no_options.txt")
output = processFile.read()
extracted = float(output[15:20])
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0=" + str(extracted) + " is outside of desired range (0.86, 0.90)"
  print "check /tmp/PERCOLATOR_no_options.txt for details" 
  success = False

# the first line of the stdout (the one after the line beginning with "PSMId")
# is extracted and if the value in the 4th column (posterior_error_prob) is 
# greater than 10e-10 an error is reported
print "(*): checking posterior_error_prob in Percolator's output..."
processFile = open("/tmp/PERCOLATOR_no_options.txt")
output = ""
line = processFile.readline()
finished = False
while (not finished): # reading line by line, looking for "PSMId"
  if line[0:5] == "PSMId":
    output = processFile.readline() # right line extracted
    finished = True
  else:
    line = processFile.readline()
tabsNumber = 0
extracted = ""
for i in range(0, len(output)): # reading char by char, looking for 4th column
  if output[i] == "\t":
    tabsNumber = tabsNumber + 1
  if tabsNumber == 3:
    extracted += output[i] # right column extracted
extracted = float(extracted)
threshold = pow(10,-10)
if extracted > threshold:
  print "...TEST FAILED: posterior_error_prob=" + str(extracted) + " is too high"
  print "check /tmp/PERCOLATOR_no_options.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)
