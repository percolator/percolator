# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of percolator
# Parameters: none

import os
import sys
import tempfile

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
pathToOutputData = "@pathToOutputData@"
xmlSupport = @xmlSupport@

class Tester:
  success = True
  failures = 0
  def doTest(self, result):
    if result:
      print("...succeeded")
    else:
      self.failures += 1
    self.success = result and self.success

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

def canPercRunThisXml(testName,flags,testFile):
  return canPercRunThis(testName,flags,testFile,"-k")

def canPercRunThisTab(testName,flags,testFile):
  return canPercRunThis(testName,flags,testFile,"",False)

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


T = Tester()

if xmlSupport:
  print("- PERCOLATOR PIN XML FORMAT")

  # running percolator to calculate psm probabilities
  print("(*) running percolator to calculate psm probabilities...")
  T.doTest(canPercRunThisXml("psms","-U","percolator/pin/pin.xml"))

  # running percolator to calculate peptide probabilities
  print("(*) running percolator to calculate peptide probabilities...")
  T.doTest(canPercRunThisXml("peptides","","percolator/pin/pin.xml"))

  # running percolator to calculate protein probabilities
  print("(*) running percolator to calculate protein probabilities...")
  T.doTest(canPercRunThisXml("proteins","-A","percolator/pin/pin.xml"))

  # running percolator to calculate psm probabilities given input file generated from sqt2pin with ptms (broken xml)
  #print("(*) running percolator on input with ptms...")
  #success = canPercRunThis("psm_ptms","-U","percolator/pin_ptms/pin_ptms.xml","-k")

  # running percolator on pin.xml with -D 4 option (description of correct features)
  print("(*) running percolator with description of correct features option...")
  T.doTest(canPercRunThisXml("D4on","-D 4 -U","percolator/pin/pin.xml"))

  # running percolator with option to generate tab-delimited input
  print("(*) running percolator to generate tab-delimited input...")
  tabData=os.path.join(pathToOutputData, "percolator/tab/percolatorTab ")
  T.doTest(canPercRunThis("tab_generate","-U -J " + tabData,"percolator/pin/pin.xml","-k",False))

print("- PERCOLATOR TAB FORMAT")

# running percolator to calculate psm probabilities
print("(*) running percolator to calculate psm probabilities...")
T.doTest(canPercRunThisTab("tab_psms","-U","percolator/tab/percolatorTab"))

# running percolator with option to process tab-delimited input
print("(*) running percolator to calculate peptide probabilities...")
T.doTest(canPercRunThisTab("tab_peptides","","percolator/tab/percolatorTab"))

# running percolator with option to process tab-delimited input
print("(*) running percolator to calculate protein probabilities...")
T.doTest(canPercRunThisTab("tab_proteins","-A","percolator/tab/percolatorTab"))

# running percolator with option to process tab-delimited input
print("(*) running percolator with description of correct features option...")
T.doTest(canPercRunThisTab("tab_D4on","-D 4 -U","percolator/tab/percolatorTabDOC"))

# if no errors were encountered, succeed
if T.success:
 print("...ALL TESTS SUCCEEDED")
 exit(0)
else:
 print("..." + T.failures + " TESTS FAILED")
 exit(1)
