# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of qvality
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
tmpDir = "@pathToData@"
success = True

print "QVALITY CORRECTNESS"

# running qvality
print "(*): running qvality..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "qvality ") + 
  os.path.join(pathToData, "qvality/target.xcorr ") + 
  os.path.join(pathToData, "qvality/null.xcorr ") + 
  "2>&1) > " + tmpDir + "/qvalityOutput.txt")

exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: qvality terminated with " + str(exitStatus) + " exit status"
  print "check qvalityOutput.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)
