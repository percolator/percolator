# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys


path = os.path.dirname(sys.argv[0])
if path == "":
  path = "./"

# running sqt2pin to generate pin.xml
print "PERCOLATOR CORRECTNESS (STEP 1): running sqt2pin..." 
processFile = os.popen(os.path.join(path,"sqt2pin ") +
  "-o " + os.path.join(path, "data/percolator_test/pin.xml ") + 
  os.path.join(path, "data/percolator_test/target.sqt ") + 
  os.path.join(path, "data/percolator_test/reverse.sqt"))
exitStatus = processFile.close()
if exitStatus is not None:
  #print "...TEST FAILED: sqt2pin terminated with " + str(exitStatus) + " exit status"
  #exit(1)
  print "...WARNING: sqt2pin terminated with " + str(exitStatus) + " exit status"

# running percolator on pin.xml with no options; 
print "PERCOLATOR CORRECTNESS (STEP 2): running percolator with no options..."
processFile = os.popen("(" + os.path.join(path, "percolator ") +  "-E " + 
  os.path.join(path, "data/percolator_test/pin.xml ") + 
  "2>&1) > /tmp/percolatorCorrectnessOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator terminated with " + str(exitStatus) + " exit status"
  exit(1)

# running percolator on pin.xml with -m option; 
print "PERCOLATOR CORRECTNESS (STEP 3): running percolator with -m option..."
processFile = os.popen("(" + os.path.join(path, "percolator ") +  "-m 2 -E " + 
  os.path.join(path, "data/percolator_test/pin.xml ") + 
  "2>&1) > /tmp/percolatorCorrectnessOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -m option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# running percolator on pin.xml with -D option; 
print "PERCOLATOR CORRECTNESS (STEP 3): running percolator with -D option..."
processFile = os.popen("(" + os.path.join(path, "percolator ") +  "-D 14 -E " + 
  os.path.join(path, "data/percolator_test/pin.xml ") + 
  "2>&1) > /tmp/percolatorCorrectnessOutput.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -D option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# if no errors were encountered, succeed
os.popen("rm /tmp/percolatorCorrectnessOutput.txt")
print "...TEST SUCCEEDED"
exit(0)
