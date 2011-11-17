# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of the percolator converters
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
success = True

print "CONVERTERS CORRECTNESS"

# SQT2PIN
# running sqt2pin with no options to generate pin.xml
# -o /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml /scratch/temp/bin/data/percolator_test/sqt2pin/target.sqt /scratch/temp/bin/data/percolator_test/sqt2pin/reverse.sqt
print "(*): running sqt2pin with no options to generate pin input..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "sqt2pin ") +
  os.path.join(pathToData, "converters/sqt2pin/target.sqt ") +
  os.path.join(pathToData, "converters/sqt2pin/reverse.sqt ") +
  "2>&1) > /tmp/SQT2PIN_no_options.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: sqt2pin with no options terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/SQT2PIN_no_options.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)

