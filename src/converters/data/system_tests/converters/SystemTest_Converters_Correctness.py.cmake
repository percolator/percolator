# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of the percolator converters
# Parameters: none

import os
import sys
import csv
import subprocess

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
pathToOutputData = "@pathToOutputData@"

class Tester:
  failures = 0
  def doTest(self, result):
    if result:
      print("...succeeded")
    else:
      self.failures += 1


def runCmd(cmd):
  result = subprocess.run(cmd, shell=True)
  return result.returncode == 0

# validating output against schema
def validate(pinXmlFile, expectedResult = True):
  print("(*): validating xml input...")
  pinSchema = doubleQuote(os.path.join(pathToData,"../../xml/percolator_in.xsd"))
  ran_successfully = runCmd("xmllint --noout --schema " + pinSchema + " " + pinXmlFile)
  if ran_successfully != expectedResult:
    print("...TEST FAILED: produced an invalid percolator xml input file")
    print("check "+pinXmlFile+" for details")
    return False
  return True

def checkNumLines(pinTabFile, expectedNumLines, expectedResult = True):
  print("(*): checking number of lines in %s..." % pinTabFile)
  numLines = sum(1 for line in open(pinTabFile, 'r'))
  result = (numLines == expectedNumLines)
  if result != expectedResult:
    print("...TEST FAILED: number of lines not as expected: %d vs. %d" % (numLines, expectedNumLines))
    return False
  return True

def checkNumTargetsAndDecoys(pinTabFile, expectedTargets, expectedDecoys, expectedResult = True):
  print("(*): checking number of target and decoy PSMs in %s..." % pinTabFile)
  reader = csv.reader(open(pinTabFile, 'r'), delimiter = '\t')
  numTargets, numDecoys = 0, 0
  for row in reader:
    if row[0] != "SpecId" and row[0] != "DefaultDirection" and len(row) >= 2:
      if int(row[1]) == 1:
        numTargets += 1
      elif int(row[1]) == -1:
        numDecoys += 1
  result = (numTargets == expectedTargets) and (numDecoys == expectedDecoys)
  if result != expectedResult:
    print("...TEST FAILED: number of lines not as expected: (%d, %d) vs. (%d, %d)" % (numTargets, numDecoys, expectedTargets, expectedDecoys))
    return False
  return True
  
def hasColumnHeader(pinTabFile, colName, expectedResult = True):
  print("(*): checking if header %s is present in %s..." % (colName, pinTabFile))
  with open(pinTabFile, 'r') as f:
    headers = f.readline()
    result = (colName in headers.split("\t"))
    if result != expectedResult:
      print("...TEST FAILED: column name %s not present." % colName)
      return False
  return True
  
# puts double quotes around the input string, needed for windows shell
def doubleQuote(path):
  return ''.join(['"',path,'"'])

def runTest(binary, testName, extraOptions = "", expectedResult = True):
  if binary == "sqt2pin":
    ext = "sqt"
  elif binary == "msgf2pin":
    ext = "mzid"
  elif binary == "tandem2pin":
    ext = "t.xml"
  else:
    print("Unknown binary %s" % binary)
    return False
  
  if expectedResult == False:
    ext = "bogus"
  
  print("(*): running %s with %s..." % (binary, testName))
  
  if testName == "metafile":
    with open(os.path.join(pathToOutputData, "target_metafile.%s.txt" % (binary)), 'w') as f:
      f.write(os.path.join(pathToData, "converters/%s/target.%s" % (binary, ext)))
    
    with open(os.path.join(pathToOutputData, "decoy_metafile.%s.txt" % (binary)), 'w') as f:
      f.write(os.path.join(pathToData, "converters/%s/decoy.%s" % (binary, ext)))
      
    cmd = ' '.join([doubleQuote(os.path.join(pathToBinaries, binary)),
      doubleQuote(os.path.join(pathToOutputData, "target_metafile.%s.txt" % (binary))),
      doubleQuote(os.path.join(pathToOutputData, "decoy_metafile.%s.txt" % (binary))),
      extraOptions,
      "2>&1 >", 
      doubleQuote(os.path.join(pathToOutputData, "%s_%s.txt" % (binary,testName)))])
  else:
    cmd = ' '.join([doubleQuote(os.path.join(pathToBinaries, binary)),
      doubleQuote(os.path.join(pathToData, "converters/%s/target.%s" % (binary, ext))),
      doubleQuote(os.path.join(pathToData, "converters/%s/decoy.%s" % (binary, ext))),
      extraOptions,
      "2>&1 >", 
      doubleQuote(os.path.join(pathToOutputData, "%s_%s.txt" % (binary,testName)))])
    
  ran_successfully = runCmd(cmd)
  if ran_successfully != expectedResult:
    print(cmd)
    print("...TEST FAILED: %s with %s terminated with %s exit status" % (binary, testName, str(exitStatus)) )
    return False
  
  testName += "_combined"
  extraOptions += " -P decoy "
  if "metafile" in testName:
    with open(os.path.join(pathToOutputData, "combined_metafile.%s.txt" % (binary)), 'w') as f:
      f.write(os.path.join(pathToData, "converters/%s/combined.%s" % (binary, ext)))
    
    cmd = ' '.join([doubleQuote(os.path.join(pathToBinaries, binary)),
      doubleQuote(os.path.join(pathToOutputData, "combined_metafile.%s.txt" % (binary))),
      extraOptions,
      "2>&1 >", 
      doubleQuote(os.path.join(pathToOutputData, "%s_%s.txt" % (binary,testName)))])
  else:
    cmd = ' '.join([doubleQuote(os.path.join(pathToBinaries, binary)),
      doubleQuote(os.path.join(pathToData, "converters/%s/combined.%s" % (binary, ext))),
      extraOptions,
      "2>&1 >", 
      doubleQuote(os.path.join(pathToOutputData, "%s_%s.txt" % (binary,testName)))])
  ran_successfully = runCmd(cmd)
  if ran_successfully != expectedResult:
    print(cmd)
    print("...TEST FAILED: %s with %s terminated with %s exit status" % (binary, testName, str(exitStatus)) )
    return False
    
  return True
  
print("CONVERTERS CORRECTNESS")

T = Tester()

ms2FileOption = "-2 " + doubleQuote(os.path.join(pathToData, "converters/small.ms2"))
xmlOutputOption = "-k " + doubleQuote(os.path.join(pathToOutputData, "%s.pin.xml"))

binaries = ["sqt2pin", "msgf2pin", "tandem2pin"]
numTargetsSeparate = [273, 203, 254]
numDecoysSeparate = [272, 215, 254]
numTargetsCombined = [169, 136, 161]
numDecoysCombined = [106, 77, 99]

for binary, nt1, nd1, nt2, nd2 in zip(binaries, numTargetsSeparate, numDecoysSeparate, numTargetsCombined, numDecoysCombined):
  print("")
  
  # verify that the tests can fail
  T.doTest(runTest(binary, "no_input_files", expectedResult = False))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "no_input_files"))
  #T.doTest(checkNumLines(pinTabFile, n+2, expectedResult = False)) # add 2 for header and default direction lines
  T.doTest(hasColumnHeader(pinTabFile, "RT", expectedResult = False))
  T.doTest(validate(pinTabFile, expectedResult = False))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt1, nd1, expectedResult = False))
  
  # run without any options
  T.doTest(runTest(binary, "no_options"))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "no_options"))
  #T.doTest(checkNumLines(pinTabFile, n+2)) # add 2 for header and default direction lines
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt1, nd1))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "no_options_combined"))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt2, nd2))
  
  # run without any options using metafile
  T.doTest(runTest(binary, "metafile"))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "metafile"))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt1, nd1))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "metafile_combined"))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt2, nd2))
  
  # run with option to add retention times as an extra column for "DOC" option in percolator
  T.doTest(runTest(binary, "RT", ms2FileOption))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "RT"))
  #T.doTest(checkNumLines(pinTabFile, n+2)) # add 2 for header and default direction lines
  T.doTest(hasColumnHeader(pinTabFile, "RT"))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt1, nd1))
  pinTabFile = os.path.join(pathToOutputData, "%s_%s.txt" % (binary, "RT_combined"))
  T.doTest(checkNumTargetsAndDecoys(pinTabFile, nt2, nd2))
  
  # run with option to output percolator input xml
  T.doTest(runTest(binary, "XML", xmlOutputOption % binary))
  T.doTest(validate(xmlOutputOption[3:] % binary))
  
  print("")

# if no errors were encountered, succeed
if T.failures == 0:
  print("...ALL TESTS SUCCEEDED")
  exit(0)
else:
  print("..." + str(T.failures) + " TESTS FAILED")
  exit(1)
