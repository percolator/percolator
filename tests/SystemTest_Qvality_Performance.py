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
print extracted
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0 is outside of desired range (0.86, 0.90)"
  exit(1)
