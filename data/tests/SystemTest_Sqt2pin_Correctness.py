# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys


path = os.path.join(os.path.dirname(sys.argv[0]), "../../")

print "SQT2PIN CORRECTNESS"

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

# uncomment to delete tests output files
#os.popen("rm /tmp/PERCOLATOR_sqt2pin.txt")


# if no errors were encountered, succeed
print "...TEST SUCCEEDED"
exit(0)
