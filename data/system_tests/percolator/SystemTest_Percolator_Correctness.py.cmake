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
  failures = 0
  def doTest(self, result):
    if result:
      print("...succeeded")
    else:
      self.failures += 1

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
  cmd = ' '.join([percExe, testFileFlag, readPath, '-S 2 -X', xmlOutput, flags, '>', txtOutput,'2>&1'])
  processFile = os.popen(cmd)
  exitStatus = processFile.close()
  if exitStatus is not None:
    print(cmd)
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

  print("(*) running percolator to calculate psm probabilities...")
  T.doTest(canPercRunThisXml("psms","-y -U","percolator/pin/pin.xml"))

  print("(*) running percolator to calculate peptide probabilities...")
  T.doTest(canPercRunThisXml("peptides","-y","percolator/pin/pin.xml"))

  print("(*) running percolator to calculate protein probabilities with picked-protein...")
  T.doTest(canPercRunThisXml("proteins","-f auto -P decoy_","percolator/pin/pin.xml"))
  
  print("(*) running percolator to calculate protein probabilities with fido...")
  T.doTest(canPercRunThisXml("proteins-fido","-A -P decoy_","percolator/pin/pin.xml"))

  print("(*) running percolator with subset training option...")
  T.doTest(canPercRunThisXml("subset_training","-y -N 1000 -U","percolator/pin/pin.xml"))
  
  print("(*) running percolator to generate tab-delimited input...")
  tabData=os.path.join(pathToOutputData, "percolatorTab ")
  T.doTest(canPercRunThis("tab_generate","-y -U -J " + tabData,"percolator/pin/pin.xml","-k",False))

# running percolator with option to process tab-delimited input
print("- PERCOLATOR TAB FORMAT")

print("(*) running percolator to calculate psm probabilities...")
T.doTest(canPercRunThisTab("tab_psms","-y -U","percolator/tab/percolatorTab"))

print("(*) running percolator to calculate peptide probabilities...")
T.doTest(canPercRunThisTab("tab_peptides","-y","percolator/tab/percolatorTab"))

print("(*) running percolator to calculate protein probabilities with picked-protein...")
T.doTest(canPercRunThisTab("tab_proteins","-f auto -P decoy_","percolator/tab/percolatorTab"))

print("(*) running percolator to calculate protein probabilities with fido...")
T.doTest(canPercRunThisTab("tab_proteins-fido","-A -P decoy_","percolator/tab/percolatorTab"))

print("(*) running percolator with subset training option...")
T.doTest(canPercRunThisTab("tab_subset_training","-y -N 1000 -U","percolator/tab/percolatorTab"))

print("(*) running percolator with static model option...")
weights = os.path.join(tempfile.gettempdir(), "test_weights.txt")
T.doTest(canPercRunThisTab("save_weights", "-w %s" % weights, "percolator/tab/percolatorTab"))
T.doTest(canPercRunThisTab("static_model", "--static -W %s" % weights, "percolator/tab/percolatorTab"))
os.remove(weights)

# if no errors were encountered, succeed
if T.failures == 0:
  print("...ALL TESTS SUCCEEDED")
  exit(0)
else:
  print("..." + str(T.failures) + " TESTS FAILED")
  exit(1)
