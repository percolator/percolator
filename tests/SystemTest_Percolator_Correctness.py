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
  print "...TEST FAILED: sqt2pin terminated with " + str(exitStatus) + " exit status"
  exit(1)

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

# running percolator on pin.xml with -D 4 option; 
# the output line containing "New pi_0 estimate" is extracted in percolator's
# output with and without -D 4 option: if the estimate with the option is lower
# than without it, an error is reported
print "PERCOLATOR CORRECTNESS (STEP 4): running percolator with -D 4 option..."
processFile = os.popen("(" + os.path.join(path, "percolator ") +  "-D 4 -E " + 
  os.path.join(path, "data/percolator_test/sqt2pin_retention_Output.xml ") + 
  "2>&1) > /tmp/percolatorOutput_D4on.txt")
exitStatus = processFile.close()
processFile = os.popen("(" + os.path.join(path, "percolator ") +  "-E " + 
  os.path.join(path, "data/percolator_test/sqt2pin_retention_Output.xml ") + 
  "2>&1) > /tmp/percolatorOutput_D4off.txt")
exitStatus = processFile.close()
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/percolatorOutput_D4on.txt")
output = processFile.read()
extracted_D4on = int(output[39:42])
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/percolatorOutput_D4off.txt")
output = processFile.read()
extracted_D4off = int(output[39:42])
if extracted_D4on < extracted_D4off:
  print "...TEST FAILED: percolator with -D 4 option performed worse than without it"
  exit(1)

# if no errors were encountered, succeed
os.popen("rm /tmp/percolatorCorrectnessOutput.txt")
os.popen("rm /tmp/percolatorOutput_D4on.txt")
os.popen("rm /tmp/percolatorOutput_D4off.txt")
print "...TEST SUCCEEDED"
exit(0)
