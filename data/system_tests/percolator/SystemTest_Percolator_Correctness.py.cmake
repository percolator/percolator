# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys
import tempfile

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
pathToOutputData = "@pathToOutputData@"
success = True

print("PERCOLATOR CORRECTNESS")

# validating output against schema
def validate(what,file):
  success = True
  print("(*) validating "+what+" output...")
  outSchema = doubleQuote(os.path.join(pathToData,"../src/xml/percolator_out.xsd"))
  processFile = os.popen("xmllint --noout --schema " + outSchema + " " + file)
  exitStatus = processFile.close()
  if exitStatus is not None:
    print("...TEST FAILED: percolator (on "+what+" probabilities) produced an invalid output")
    print("check "+file+" for details")
    success = False
  return success

def canPercRunThis(testName,flags,testFile,testFileFlag="",checkValidXml=True):
  success = True
  outputPath = os.path.join(pathToOutputData,"PERCOLATOR_"+testName)
  xmlOutput = doubleQuote(outputPath + ".pout.xml")
  txtOutput = doubleQuote(outputPath + ".txt")
  readPath = doubleQuote(os.path.join(pathToData, testFile))
  percExe = doubleQuote(os.path.join(pathToBinaries, "percolator"))
  processFile = os.popen(' '.join([percExe, testFileFlag, readPath, '-X', xmlOutput, flags, '>', txtOutput,'2>&1']))
  exitStatus = processFile.close()
  if exitStatus is not None:
    print(' '.join([percExe, testFileFlag, readPath, '-X', xmlOutput, flags, '2>&1 >', txtOutput]))
    print("...TEST FAILED: percolator ("+testName+") terminated with " + os.strerror(exitStatus) + " exit status")
    print("check "+ txtOutput +" for details") 
    success = False
  if checkValidXml:
    success = success and validate(testName,xmlOutput)
  return success

# puts double quotes around the input string, needed for windows shell
def doubleQuote(path):
  return ''.join(['"',path,'"'])

# running percolator to calculate psm probabilities
print("(*) running percolator to calculate psm probabilities...")
success = canPercRunThis("psms","-U","percolator/pin/pin.xml","-k") and success

# running percolator to calculate peptide probabilities
print("(*) running percolator to calculate peptide probabilities...")
success = canPercRunThis("peptides","","percolator/pin/pin.xml","-k") and success

# running percolator to calculate protein probabilities
print("(*) running percolator to calculate protein probabilities...")
success = canPercRunThis("proteins","-A","percolator/pin/pin.xml","-k") and success

# running percolator to calculate psm probabilities given input file generated from sqt2pin with ptms
#print("(*) running percolator on input with ptms...")
#success = canPercRunThis("psm_ptms","-U","percolator/pin_ptms/pin_ptms.xml","-k")

# running percolator with option to generate tab-delimited input
print("(*) running percolator to generate tab-delimited input...")
tabData=os.path.join(pathToData, "percolator/tab/percolatorTab ")
success = canPercRunThis("tab_generate","-U -J " + tabData,"percolator/pin/pin.xml","-k",False) and success

# running percolator with option to process tab-delimited input
print("(*) running percolator on tab-delimited input...")
success = canPercRunThis("tab","-U","percolator/tab/percolatorTab","",False) and success

# running percolator on pin.xml with -D 4 option (description of correct features)
print("(*) running percolator with description of correct features option...")
success = canPercRunThis("D4on","-D 4 -U","percolator/pin/pin.xml","-k") and success

# if no errors were encountered, succeed
if success == True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)
