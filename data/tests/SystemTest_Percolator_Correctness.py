# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys


path = os.path.join(os.path.dirname(sys.argv[0]), "../../")

print "PERCOLATOR CORRECTNESS"

# PERCOLATOR_sqt2pin
# running percolator on pin.xml with no options; 
# -E /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml -X /scratch/temp/bin/data/percolator_test/sqt2pin/pout.xml
print "(STEP 1): running percolator on pin input with no options..."
processFile = os.popen("(" +"percolator " + "-E " + 
  os.path.join(path, "data/percolator_test/sqt2pin/pin.xml ") + "-X " +
  os.path.join(path, "data/percolator_test/sqt2pin/pout.xml ") +
  "2>&1) > /tmp/PERCOLATOR_sqt2pin.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator (no options) terminated with " + str(exitStatus) + " exit status"
  exit(1)

processFile = os.popen("xmllint --noout --schema http://per-colator.com/xml/xml-1-1/percolator_out.xsd " +
  os.path.join(path, "data/percolator_test/sqt2pin/pout.xml"))
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator (no options) produced an invalid/ill-formed output: " + os.path.join(path, "data/percolator_test/sqt2pin/pout.xml")
  exit(1)

# PERCOLATOR_m
# running percolator on pin.xml with -m option;
# -m 2 -E /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml -X /scratch/temp/bin/data/percolator_test/m/pout.xml
print "(STEP 2): running percolator on pin input with -m option..."
processFile = os.popen("(" + "percolator " + "-m 2 -E " + 
  os.path.join(path, "data/percolator_test/sqt2pin/pin.xml ") + "-X " +
  os.path.join(path, "data/percolator_test/m/pout.xml ") +
  "2>&1) > /tmp/PERCOLATOR_m.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -m option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# PERCOLATOR_tab-delimited_generate
# running percolator with option to generate tab-delimited input
# -E /scratch/temp/bin/data/percolator_test/sqt2pin/pin.xml -J /scratch/temp/bin/data/percolator_test/tab-delimited/percolatorTab
print "(STEP 3): running percolator to generate tab-delimited input..."
processFile = os.popen("(" + "percolator " +  "-E " + 
  os.path.join(path, "data/percolator_test/sqt2pin/pin.xml ") +   "-J " +
  os.path.join(path, "data/percolator_test/tab-delimited/percolatorTab ") + 
  "2>&1) > /tmp/PERCOLATOR_tab-delimited_generate.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -J option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# PERCOLATOR_tab-delimited
# running percolator with option to process tab-delimited input
# -j /scratch/temp/bin/data/percolator_test/tab-delimited/percolatorTab -X /scratch/temp/bin/data/percolator_test/tab-delimited/pout.xml
print "(STEP 4): running percolator on tab-delimited input..."
processFile = os.popen("(" + "percolator " +  "-j " + 
  os.path.join(path, "data/percolator_test/tab-delimited/percolatorTab ") + "-X " +
  os.path.join(path, "data/percolator_test/tab-delimited/pout.xml ") + 
  "2>&1) > /tmp/PERCOLATOR_tab-delimited.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -j option terminated with " + str(exitStatus) + " exit status"
  exit(1)

# PERCOLATOR_retention_truncated
# running percolator on pin.xml with -D 4 option; 
# the output line containing "New pi_0 estimate" is extracted in percolator's
# output with and without -D 4 option: if the estimate with the option is lower
# than without it, an error is reported
# -E /scratch/temp/bin/data/percolator_test/retention_truncated/pin.xml -X /scratch/temp/bin/data/percolator_test/retention_truncated/pout_D4off.xml
# -D 4 -E /scratch/temp/bin/data/percolator_test/retention_truncated/pin.xml -X /scratch/temp/bin/data/percolator_test/retention_truncated/pout_D4on.xml
print "(STEP 5): running percolator with -D 4 option (retention times)..."
processFile = os.popen("(" + "percolator " +  "-D 4 -E " + 
  os.path.join(path, "data/percolator_test/retention_truncated/pin.xml ") +  "-X " +
  os.path.join(path, "data/percolator_test/retention_truncated/pout_D4on.xml ") +
  "2>&1) > /tmp/PERCOLATOR_retention_truncated_D4on.txt")
exitStatus = processFile.close()
processFile = os.popen("(" + "percolator " + "-E " + 
  os.path.join(path, "data/percolator_test/retention_truncated/pin.xml ") + "-X " +
  os.path.join(path, "data/percolator_test/retention_truncated/pout_D4off.xml ") + 
  "2>&1) > /tmp/PERCOLATOR_retention_truncated_D4off.txt")
exitStatus = processFile.close()
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/PERCOLATOR_retention_truncated_D4on.txt")
output = processFile.read()
extracted_D4on = int(output[39:40])
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/PERCOLATOR_retention_truncated_D4off.txt")
output = processFile.read()
extracted_D4off = int(output[39:40])
if extracted_D4on < extracted_D4off:
  print "...TEST FAILED: percolator with -D 4 option performed worse than without it"
  exit(1)

# uncomment to delete tests output files
#os.popen("rm /tmp/PERCOLATOR_m.txt")
#os.popen("rm /tmp/PERCOLATOR_tab-delimited_generate.txt")
#os.popen("rm /tmp/PERCOLATOR_tab-delimited.txt")
#os.popen("rm /tmp/PERCOLATOR_retention_truncated_D4on.txt")
#os.popen("rm /tmp/PERCOLATOR_retention_truncated_D4off.txt")

# if no errors were encountered, succeed
print "...TEST SUCCEEDED"
exit(0)
