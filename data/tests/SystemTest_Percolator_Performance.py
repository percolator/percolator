# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys


path = os.path.join(os.path.dirname(sys.argv[0]), "../../")

print "PERCOLATOR PERFORMANCE"

# SQT2PIN
# running sqt2pin to generate pin.xml
# -o /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml /scratch/temp/bin/data/percolator_test/sqt2pin/target.sqt /scratch/temp/bin/data/percolator_test/sqt2pin/reverse.sqt
print "(STEP 1): running sqt2pin to generate pin input..." 
processFile = os.popen("sqt2pin " +
  "-o " + os.path.join(path, "data/percolator_test/sqt2pin/pin.xml ") + 
  os.path.join(path, "data/percolator_test/sqt2pin/target.sqt ") + 
  os.path.join(path, "data/percolator_test/sqt2pin/reverse.sqt"))
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED:\n" + command + "\nterminated with " + str(exitStatus) + " exit status"
  exit(1)

# PERCOLATOR_sqt2pin
# running percolator on pin.xml with no options; 
# -E /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml -X /scratch/temp/bin/data/percolator_test/sqt2pin/pout.xml
print "(STEP 2): running percolator on pin input to generate pout..."
processFile = os.popen("(" + "percolator " + "-E " + 
  os.path.join(path, "data/percolator_test/sqt2pin/pin.xml ") + "-X " +
  os.path.join(path, "data/percolator_test/sqt2pin/pout.xml ") +
  "2>&1) > /tmp/PERCOLATOR_sqt2pin.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator (no options) terminated with " + str(exitStatus) + " exit status"
  exit(1)

# the output line containing "New pi_0 estimate" is extracted and if its value is 
# outside of 622+/-5% an error is reported
print "(STEP 3): checking new pi_0 estimate in pout..."
processFile = os.popen("grep \"New pi_0 estimate\" " + 
  "/tmp/PERCOLATOR_sqt2pin.txt")
output = processFile.read()
extracted = float(output[39:42])
if extracted < 590.9 or extracted > 653.1: 
  print "...TEST FAILED: new pi_0 estimate on merged list outside of desired range (590.9, 653.1)"
  exit(1)

# the output line containing "Selecting pi_0" is extracted and if its value is 
# outside of (0.86, 0.90) an error is reported
print "(STEP 4): checking selected new pi_0 estimate in pout..."
processFile = os.popen("grep \"Selecting pi_0\" " +
  "/tmp/PERCOLATOR_sqt2pin.txt")
output = processFile.read()
extracted = float(output[15:20])
if extracted < 0.86 or extracted > 0.90:
  print "...TEST FAILED: selected pi_0 is outside of desired range (0.86, 0.90)"
  exit(1)

# the first line of the stdout (the one after the line beginning with "PSMId")
# is extracted and if the value in the 4th column (posterior_error_prob) is 
# greater than 10e-10 an error is reported
print "(STEP 5): checking posterior_error_prob in pout..."
processFile = open("/tmp/PERCOLATOR_sqt2pin.txt")
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

# uncomment to delete tests output files
#os.popen("rm /tmp/PERCOLATOR_sqt2pin.txt")

# if no errors were encountered, succeed
print "...TEST SUCCEEDED"
exit(0)
