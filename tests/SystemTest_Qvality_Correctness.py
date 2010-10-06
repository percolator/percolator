# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of qvality
# Parameters: none

import os
import sys


path = os.path.dirname(sys.argv[0])
if path == "":
  path = "./"

print "QVALITY CORRECTNESS"

# running qvality
print "(STEP 1): running qvality..."
processFile = os.popen("(" + os.path.join(path, "qvality ") + 
  os.path.join(path, "data/qvality_test/target.xcorr ") + 
  os.path.join(path, "data/qvality_test/null.xcorr ") + 
  "2>&1) > /tmp/qvalityCorrectnessOutput.txt")

exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: qvality terminated with " + str(exitStatus) + " exit status"
  exit(1)

# if no errors were encountered, succeed
#os.popen("rm /tmp/qvalityCorrectnessOutput.txt")
print "...TEST SUCCEEDED"
exit(0)
