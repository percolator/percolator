# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
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

# running percolator on pin.xml with no options; 
print "PERCOLATOR PERFORMANCE (STEP 2): running percolator with no options..."
processFile = os.popen("(" + path + "/percolator " + path + 
  "-E /data/percolator_test/pin.xml 2>&1) > /tmp/percolatorPerformanceOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator terminated with " + str(exitStatus) + " exit status"
  exit(1)

# running percolator on pin.xml with -m option; 
print "PERCOLATOR PERFORMANCE (STEP 3): running percolator with -m option..."
processFile = os.popen("(" + path + "/percolator " + path + 
  "-m 2 -E /data/percolator_test/pin.xml 2>&1) > /tmp/percolatorPerformanceOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -m option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# running percolator on pin.xml with -D option; 
print "PERCOLATOR PERFORMANCE (STEP 3): running percolator with -D option..."
processFile = os.popen("(" + path + "/percolator " + path + 
  "-D 14 -E /data/percolator_test/pin.xml 2>&1) > /tmp/percolatorPerformanceOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -D option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# if no errors were encountered, succeed
os.popen("rm /tmp/percolatorPerformanceOutput.txt")
print "...TEST SUCCEEDED"
exit(0)
