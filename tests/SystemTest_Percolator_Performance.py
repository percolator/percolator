# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys


path = os.path.dirname(sys.argv[0])

# running sqt2pin to generate pin.xml
print "PERCOLATOR PERFORMANCE (STEP 1): running sqt2pin..." 
processFile = os.popen(path + "/sqt2pin -o " + path + 
  "/data/percolator_test/pin.xml " + path + "/data/percolator_test/target.sqt " 
  + path + "/data/percolator_test/reverse.sqt")
exitStatus = processFile.close()
#if exitStatus is not None:
#  print "...TEST FAILED: sqt2pin terminated with " + str(exitStatus) + 
#  " exit status"
#  exit(1)

# running percolator on pin.xml; 
print "PERCOLATOR PERFORMANCE (STEP 2): running percolator..."
processFile = os.popen("(" + path + "/percolator " + path + 
  "/data/percolator_test/pin.xml 2>&1) > /tmp/percolatorPerformanceOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator terminated with " + str(exitStatus) + " exit status"
  exit(1)

# the output line containing "New pi_0 estimate" is extracted and if its value is 
# outside of 622+/-5% an error is reported
processFile = os.popen("grep \"New pi_0 estimate\" " + 
  "/tmp/percolatorPerformanceOutput.txt")
output = processFile.read()
extracted = float(output[39:42])
if extracted < 590.9 or extracted > 653.1: 
  print "...TEST FAILED: pi_0 estimate on merged list outside of desired range (590.9, 653.1)"
  exit(1)

# the output line containing "Selecting pi_0" is extracted and if its value is 
# outside of (0.86, 0.90) an error is reported
processFile = os.popen("grep \"Selecting pi_0\" " +
  "/tmp/percolatorPerformanceOutput.txt")
output = processFile.read()
extracted = float(output[15:20])
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0 is outside of desired range (0.86, 0.90)"
  exit(1)

# the first line of the stdout (the one after the line beginning with "PSMId")
# is extracted and if the value in the 4th column (posterior_error_prob) is 
# greater than 10e-10 an error is reported
processFile = open( "/tmp/percolatorPerformanceOutput.txt" )
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
  print "...TEST FAILED: posterior_error_prob is too high (" + str(extracted) + ")"
  exit(1)

os.popen("rm /tmp/percolatorPerformanceOutput.txt")
print "...TEST SUCCEEDED"
exit(0)
