# Mattia Tomasoni - Percolator Project
# Script that tests the performances of qvality
# Parameters: none

import os
import sys


path = os.path.dirname(sys.argv[0])


# running qvality
print "QVALITY PERFORMANCE (STEP 1): running qvality..."
processFile = os.popen("(" + path + "/qvality " + path + 
  "/data/qvality_test/target.xcorr " + path + 
  "/data/qvality_test/null.xcorr 2>&1) > /tmp/qvalityPerformanceOutput.txt")
exitStatus = processFile.close()
#if exitStatus is not None:
#  print "...TEST FAILED: qvality terminated with " + str(exitStatus) + " exit status"
#  exit(1)


# the output line containing "Selecting pi_0" is extracted and if its value is 
# outside of (0.86, 0.90) an error is reported
processFile = os.popen("grep \"Selecting pi_0\" " +
  "/tmp/qvalityPerformanceOutput.txt")
output = processFile.read()
extracted = float(output[15:20])
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0 is outside of desired range (0.86, 0.90)"
  exit(1)

# the number of lines of stdout (after the line beginning with "Score") until 
# q-value < 0.01 are counted and an error is reported if their number is greater
# than 755+/-5%
processFile = open("/tmp/qvalityPerformanceOutput.txt")
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
  print "...TEST FAILED: value outside of desired range (717, 793)"
  exit(1)

# if no errors were encountered, succeed
#os.popen("rm /tmp/percolatorPerformanceOutput.txt")
print "...TEST SUCCEEDED"
exit(0)
