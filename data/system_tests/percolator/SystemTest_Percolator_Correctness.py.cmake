# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
success = True

print "PERCOLATOR CORRECTNESS"

# PERCOLATOR
# running percolator with no options
# /scratch/temp/bin/data/percolator/pin/pin.xml -X /scratch/temp/bin/data/percolator/tab-delimited/pout.xml -U
print "(*) running percolator with no options..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/pin/pin.xml ") + "-X " +
  os.path.join(pathToData, "percolator/pin/pout.xml ") + "-U " +
  "2>&1) > /tmp/PERCOLATOR_no_options.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with no options terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_no_options.txt for details" 
  success = False

# PERCOLATOR_tab-delimited_generate
# running percolator with option to generate tab-delimited input
# /scratch/temp/bin/data/percolator/sqt2pin/pin.xml -J /scratch/temp/bin/data/percolator/tab/percolatorTab -U
print "(*) running percolator to generate tab-delimited input..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") + 
  os.path.join(pathToData, "percolator/pin/pin.xml ") + "-J " +
  os.path.join(pathToData, "percolator/tab/percolatorTab ") + "-U " +
  "2>&1) > /tmp/PERCOLATOR_tab_generate.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -J option terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_tab_generate.txt for details" 
  success = False

# PERCOLATOR_tab-delimited
# running percolator with option to process tab-delimited input
# -j /scratch/temp/bin/data/percolator/tab/percolatorTab -X /scratch/temp/bin/data/percolator/tab/pout.xml -U
print "(*) running percolator on tab-delimited input..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +  "-j " + 
  os.path.join(pathToData, "percolator/tab/percolatorTab ") + "-X " +
  os.path.join(pathToData, "percolator/tab/pout.xml ") + "-U " +
  "2>&1) > /tmp/PERCOLATOR_tab.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -j option terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_tab.txt for details" 
  success = False

# PERCOLATOR_rt
# running percolator on pin.xml with -D 4 option; 
# the output line containing "New pi_0 estimate" (estimated number of psms at
# q-value 0.01 threshold) is extracted from percolator's output with and 
# without -D 4 option: if the estimate with the option is lower than without 
# it, an error is reported.
# /scratch/temp/bin/data/percolator/rt/pin.xml -X /scratch/temp/bin/data/percolator/rt/pout_D4off.xml -D 4 -U
# /scratch/temp/bin/data/percolator/rt/pin.xml -X /scratch/temp/bin/data/percolator/rt/pout_D4on.xml -U
print "(*) running percolator with -D 4 option (retention times)..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/rt/pin.xml ") +  "-X " +
  os.path.join(pathToData, "percolator/rt/pout_D4on.xml ") + "-D 4 -U " + 
  "2>&1) > /tmp/PERCOLATOR_rt_D4on.txt")
exitStatus = processFile.close()
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/rt/pin.xml ") + "-X " +
  os.path.join(pathToData, "percolator/rt/pout_D4off.xml ") + "-U " + 
  "2>&1) > /tmp/PERCOLATOR_rt_D4off.txt")
exitStatus = processFile.close()
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/PERCOLATOR_rt_D4on.txt")
output = processFile.read()
extracted_D4on = int(output[39:40])
processFile = os.popen("grep \"New pi_0\" " +
  "/tmp/PERCOLATOR_rt_D4off.txt")
output = processFile.read()
extracted_D4off = int(output[39:40])
if extracted_D4on < extracted_D4off:
  print "...TEST FAILED: percolator with -D 4 option performed worse than without it"
  print "check /tmp/PERCOLATOR_rt_D4on.txt and /tmp/PERCOLATOR_rt_D4off.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)
