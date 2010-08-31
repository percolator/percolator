# @Created by L. Moruz
# 30th Aug, 2010 
# Script that tests Elude's functionality
# Parameters: none

import os
import sys

########################## TEST 1: train a model and apply it to predict rt ##########################
# this test makes sure that when using -t -e we get i) the correct rt predictions, ii) the correct performance (correlation and delta t)
# this test are performed using various formats of the input files
def testStandalone():
  print "TEST 1 (train a model)..."
  path = os.path.dirname(sys.argv[0])
  # folder where the test data is located
  testFolder = os.path.join(path, "data/elude_test/standalone/")
  # temporary folder to store the output files 
  tmpFolder = os.path.join(path, "data/elude_test/tmp/")
  # path to elude 
  eludePath = os.path.join(path, "elude")
  
  # make a temporary folder to store the output files 
  if os.path.isdir(tmpFolder):
    removeDir(tmpFolder)
  os.mkdir(tmpFolder)
 
  ############# Test 1.1 - the format of the ppetide sequence in train and test files is X.A.Y, test file does not include rt
  # run Elude  
  os.system(eludePath + " -t " + testFolder + "train.txt -e " + testFolder + "test.txt -o " + tmpFolder + "out -s " +  tmpFolder + "model -i " + tmpFolder + "insource -g " + tmpFolder + "index -y -u -x -v 5 2> " + tmpFolder + "log")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "out", "1191", 1107, "R.QLLACIVNPEIMK.T\t43.2946"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.1: Elude standalone: Incorrect predictions in the output file")

  # check that the in source fragments are correct
  if not testFileContent(tmpFolder + "insource", "2", -1, ""):
    removeDir(tmpFolder)
    failTest("Elude Test 1.1: Elude standalone: Incorrect in source fragments")

  # check the retention index 
  if not testFileContent(tmpFolder + "index", "20", 18, "V\t0.581933"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.1: Elude standalone: Incorrect retention index")

  # check the load model 
  os.system(eludePath + " -m " + tmpFolder + "model -e " + testFolder + "test.txt -o " + tmpFolder + "outLoad " + " -y -u -x -v 5 2> " + tmpFolder + "log2")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "outLoad", "1191", 1107, "R.QLLACIVNPEIMK.T\t43.2946"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.1:  Elude standalone: Incorrect predictions in the output file when loading a model")

  ############# Test 1.2 - the format of the peptide sequence in train and test files is A, test file does not include rt
  # run Elude; results are overwritten  
  os.system(eludePath + " -t " + testFolder + "train_1.txt -e " + testFolder + "test_1.txt -o " + tmpFolder + "out -s " +  tmpFolder + "model -i " + tmpFolder + "insource -g " + tmpFolder + "index -y -u -x -v 5 2> " + tmpFolder + "log")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "out", "153", 36, "GIAMPTCISIDNCICHYTPLK\t45.5099"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.2:  Elude standalone: Incorrect predictions in the output file")

  # check that the in source fragments are correct
  if not testFileContent(tmpFolder + "insource", "2", -1, ""):
    removeDir(tmpFolder)
    failTest("Elude Test 1.2:  Elude standalone: Incorrect in source fragments")

  # check the retention index 
  if not testFileContent(tmpFolder + "index", "20", 7, "H\t-0.607149"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.2: Incorrect retention index")

  # check the load model 
  os.system(eludePath + " -m " + tmpFolder + "model -e " + testFolder + "test_1.txt -o " + tmpFolder + "outLoad " + " -y -u -x -v 5 2> " + tmpFolder + "log2")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "out", "153", 36, "GIAMPTCISIDNCICHYTPLK\t45.5099"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.2: Elude standalone: Incorrect predictions in the output file when loading a model")

 ############# Test 1.3 - the format of the peptide sequence in train and test files is A.X.Y, test file includes rt
  # run Elude; results are overwritten  
  os.system(eludePath + " -t " + testFolder + "train.txt -e " + testFolder + "test_2.txt -o " + tmpFolder + "out -s " +  tmpFolder + "model -i " + tmpFolder + "insource -g " + tmpFolder + "index -y -u -x -v 5 2> " + tmpFolder + "log")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "out", "1742", 695, "K.FAILTLNSFLK.D\t53.8287\t78.2652"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.3:  Elude standalone: Incorrect predictions in the output file")

  # check that the in source fragments are correct
  if not testFileContent(tmpFolder + "insource", "3", 3, "R.EETGK.V\t26.8869\tTest"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.3:  Elude standalone: Incorrect in source fragments")

  # check the retention index 
  if not testFileContent(tmpFolder + "index", "20", 17, "T\t0.139622"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.3: Incorrect retention index")

  # check performance 
  f = open(tmpFolder + "log")
  lines = f.readlines(90)
  if (lines[86].find("0.898") == -1) or (lines[88].find("44.3") == -1): 
    removeDir(tmpFolder)
    failTest("Elude Test 1.3: Incorrect performance measures")

  # check the load model 
  os.system(eludePath + " -m " + tmpFolder + "model -e " + testFolder + "test_1.txt -o " + tmpFolder + "outLoad " + " -y -u -x -v 5 2> " + tmpFolder + "log2")

  # check the file with the output predictions
  if not testFileContent(tmpFolder + "out", "1742", 695, "K.FAILTLNSFLK.D\t53.8287\t78.2652"):
    removeDir(tmpFolder)
    failTest("Elude Test 1.3: Elude standalone: Incorrect predictions in the output file when loading a model")
  
  removeDir(tmpFolder)
  print "...TEST SUCCEEDED"
  return 0
  

########################## TEST 2: calibrate a model from the library ##########################
def testCalibration():
  print "TEST 2 (calibrate a model)..."
  path = os.path.dirname(sys.argv[0])
  # folder where the test data is located
  testFolder = os.path.join(path, "data/elude_test/calibrate_data/")
  # temporary folder to store the output files 
  tmpFolder = os.path.join(path, "data/elude_test/tmp/")
  # path to elude 
  eludePath = os.path.join(path, "elude")
  
  # make a temporary folder to store the output files 
  if os.path.isdir(tmpFolder):
    removeDir(tmpFolder)
  os.mkdir(tmpFolder)
 
  # run Elude  
  os.system(eludePath + " -t " + testFolder + "calibrate.txt" +  " -e " + testFolder + "test.txt -o " + tmpFolder + "out " + "-a -u -x -v 5 2> " + tmpFolder + "log")

  # check the performance 
  f = open(tmpFolder + "log")
  lines = f.readlines(279) 
  if (lines[275].find("0.943") == -1) or (lines[277].find("43.17") == -1):
    removeDir(tmpFolder)
    failTest("Elude Test 2: Elude calibration: Incorrect performance measures")

  # check the predictions 
  if not testFileContent(tmpFolder + "out", "1743", 1719, "K.IGDYAGIK.W\t20.8995\t22.1232"):
    removeDir(tmpFolder)
    failTest("Elude Test 2: Elude calibration: Incorrect predictions in the output file")

  removeDir(tmpFolder)
  print "...TEST SUCCEEDED"
  return 0


########################## Useful functions ##########################
# deletes the contents of a directory 
def removeDir(dirname):
  files = os.listdir(dirname)
  for f in files:
    os.remove(os.path.join(dirname, f))
  os.rmdir(dirname) 

# check if the file filename includes noLines lines and that the content if line lineNumber is lineContent
def testFileContent(filename, noLines, lineNumber, lineContent):
  # check the number of lines
  output = os.popen("wc -l " + filename + " 2>&1 | grep " + noLines,"r")
  if (output.readline() == "") or (lineNumber > noLines):
    return False 

  if lineNumber == -1:
    return True 

   # check the content of lineNumber
  handle = open(filename, "r")
  i = 1
  while i < lineNumber: 
    handle.readline()
    i = i + 1

  if handle.readline().strip() != lineContent.strip():
    return False

  return True

# fail a test with a message 
def failTest(message):
  print message 
  print "...TEST FAILED"
  exit(1)


########################## MAIN ##########################
print "\nELUDE PERFORMANCE" 

testStandalone()
testCalibration()

print "ALL TESTS COMPLETED SUCCESSFULLY\n"
exit(0)

