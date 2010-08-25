# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os

# USER MUST SET THESE VARIABLES
# target sqt file
targetSqt="/home/mattia/percolator/data/percolator_test/target.sqt" 
# the decoy sqt file
decoySqt="/home/mattia/percolator/data/percolator_test/reverse.sqt"
# xml file in which to save the output of sqt2pin 
pin="/home/mattia/percolator/data/percolator_test/pin.xml"
# percolator executable 
percolator="/home/mattia/percolatorBuild/src/percolator"
# sqt2pin executable
sqt2pin="/home/mattia/percolatorBuild/src/converters/sqt2pin" 

# running sqt2pin to generate pin.xml
print "PERCOLATOR PERFORMANCE (STEP 1): running sqt2pin..." 
os.popen(sqt2pin + " -x " + pin + " " + targetSqt + " " + decoySqt + 
  " > /dev/null")

# running percolator on pin.xml; the output line containing pi_0 is extracted
# and if its value is outside of 622+/-5% an error is reported
print "PERCOLATOR PERFORMANCE (STEP 2): running percolator..."
output = os.popen(percolator + " " + pin + " 2>&1 | grep \"New pi_0 estimate\"")
outputStr = output.read()
# extracting pi_0 value
pi_0 = int(outputStr[39:42])
# check whether pi_0 is within desired range
if pi_0 < 590.9 or pi_0 > 653.1: 
  print "...TEST FAILED"
else: 
  print "...TEST SUCCEEDED"
