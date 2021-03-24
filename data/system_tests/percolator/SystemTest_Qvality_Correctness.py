# Mattia Tomasoni - Percolator Project
# Script that tests the correctness of qvality
# Parameters: none

import os
import sys

pathToBinaries = "/usr/local/bin"
pathToData = "/home/ekvall/percolator/data"
pathToOutputData = "/home/ekvall/percolator/data"
success = True

print("QVALITY CORRECTNESS")

# puts double quotes around the input string, needed for windows shell
def doubleQuote(path):
  return ''.join(['"',path,'"'])

# running qvality
print("(*): running qvality...")
processFile = os.popen(' '.join([doubleQuote(os.path.join(pathToBinaries, "qvality")),
  doubleQuote(os.path.join(pathToData, "qvality/target.xcorr")),
  doubleQuote(os.path.join(pathToData, "qvality/null.xcorr")), 
  '>', doubleQuote(os.path.join(pathToOutputData, "qvalityOutput.txt")),'2>&1']))

exitStatus = processFile.close()
if exitStatus is not None:
  print(' '.join([doubleQuote(os.path.join(pathToBinaries, "qvality")),
    doubleQuote(os.path.join(pathToData, "qvality/target.xcorr")),
    doubleQuote(os.path.join(pathToData, "qvality/null.xcorr")), 
    '>', doubleQuote(os.path.join(pathToOutputData, "qvalityOutput.txt")),'2>&1']))
  print("...TEST FAILED: qvality terminated with " + str(exitStatus) + " exit status")
  print("check qvalityOutput.txt for details")
  success = False

# if no errors were encountered, succeed
if success == True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)
