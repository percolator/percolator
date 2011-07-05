# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
success = True

print "PERCOLATOR CORRECTNESS"

# validating output against schema
def validate(what,file):
  success = True
  print "(*) validating "+what+" output..."
  processFile = os.popen("xmllint --noout --schema " + pathToData + "/../src/xml/percolator_out.xsd " + file)
  exitStatus = processFile.close()
  if exitStatus is not None:
    print "...TEST FAILED: percolator (on "+what+" probabilities) produced an invalid output" 
    print "check "+file+" for details"
    success = False
  return success

# running percolator to calculate psm probabilities
print "(*) running percolator to calculate psm probabilities..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/pin/pin.xml ") + "-X /tmp/PERCOLATOR_psm.pout.xml " + "-U " +
  "2>&1) > /tmp/PERCOLATOR_psm.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator (psm probabilities) terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_psm.txt for details" 
  success = False
success = validate("psms","/tmp/PERCOLATOR_psm.pout.xml")

# running percolator to calculate peptide probabilities
print "(*) running percolator to calculate peptide probabilities..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/pin/pin.xml ") + "-X /tmp/PERCOLATOR_peptide.pout.xml " +
  "2>&1) > /tmp/PERCOLATOR_peptide.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator (peptide probabilities) terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_peptide.txt for details" 
  success = False
success = validate("peptides","/tmp/PERCOLATOR_peptide.pout.xml")

# running percolator with option to generate tab-delimited input
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

# running percolator with option to process tab-delimited input
print "(*) running percolator on tab-delimited input..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +  "-j " + 
  os.path.join(pathToData, "percolator/tab/percolatorTab ") + "-X /tmp/PERCOLATOR_tab.pout.xml -U " +
  "2>&1) > /tmp/PERCOLATOR_tab.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -j option terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_tab.txt for details" 
  success = False

# running percolator on pin.xml with -D 4 option (description of correct features)
print "(*) running percolator with description of correct features option..."
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/rt/pin.xml ") +  "-X /tmp/PERCOLATOR_D4on.pout.xml -D 4 -U " + 
  "2>&1) > /tmp/PERCOLATOR_D4on.txt")
exitStatus = processFile.close()
processFile = os.popen("(" + os.path.join(pathToBinaries, "percolator ") +
  os.path.join(pathToData, "percolator/rt/pin.xml ") + "-X /tmp/PERCOLATOR_D4off.pout.xml -U " + 
  "2>&1) > /tmp/PERCOLATOR_D4off.txt")
exitStatus = processFile.close()
if exitStatus is not None:
  print "...TEST FAILED: percolator with -D 4 option terminated with " + str(exitStatus) + " exit status"
  print "check /tmp/PERCOLATOR_D4on.txt and /tmp/PERCOLATOR_D4off.txt for details" 
  success = False

# if no errors were encountered, succeed
if success == True:
 print "...TEST SUCCEEDED"
 exit(0)
else:
 print "...TEST FAILED"
 exit(1)
