# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys

path = os.path.dirname(sys.argv[0])

# running sqt2pin to generate pin.xml
print "PERCOLATOR PERFORMANCE (STEP 1): running sqt2pin..." 
os.popen(path + "/sqt2pin -x " + path + "/data/percolator_test/pin.xml " + path 
  + "/data/percolator_test/target.sqt " + path + "/data/percolator_test/reverse.sqt")

# running percolator on pin.xml; the output line containing pi_0 is extracted
# and if its value is outside of 622+/-5% an error is reported
print "PERCOLATOR PERFORMANCE (STEP 2): running percolator..."
output = os.popen(path + "/percolator " + path + "/data/percolator_test/pin.xml 2>&1 | grep \"New pi_0 estimate\"")
outputStr = output.read()
# extracting pi_0 value
pi_0 = int(outputStr[39:42])
# check whether pi_0 is within desired range
if pi_0 < 590.9 or pi_0 > 653.1: 
  print "...TEST FAILED"
  exit(1)
else: 
  print "...TEST SUCCEEDED"
  exit(0)
