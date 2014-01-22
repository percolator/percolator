# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
tmpDir = "@pathToData@"
success = True

print("PERCOLATOR PERFORMANCE")

# check number of psms/peptides with q-value < 0.01 is withing 5% from expected
def checkNumberOfSignificant(what,file,expected):
  success = True
  print("(*): checking number of significant "+what+" found...")
  if what=="proteins":
    processFile = os.popen("grep \"The number of Proteins idenfified at q-value = 0.01 is :\" "+file)
  else:
    processFile = os.popen("grep \"New pi_0 estimate\" "+file)
  output = processFile.read()
  if what=="proteins":
    extracted = float(output[-3:])
  else:
    extracted = float(output[39:42])
  if extracted<expected-(5*expected/100)  or extracted>expected+(5*expected/100) : 
    print("...TEST FAILED: number of significant "+what+"=" + str(int(extracted)) + " is outside of desired range")
    print("check "+file+" for details") 
    success = False
  return success

# check value of pi_0 is within 5% expected value
def checkPi0(what,file,expected):
  success = True
  print("(*): checking pi_0 estimate for "+what+"...")
  processFile = os.popen("grep \"Selecting pi_0\" "+file)
  output = processFile.read()
  extracted = float(output[15:20])
  if extracted<expected-(5*expected/100) or extracted>expected+(5*expected/100):
    print("...TEST FAILED: "+what+" pi_0=" + str(extracted) + " is outside of desired range")
    print("check "+file+" for details") 
    success = False
  return success

# check pep within 5% expected value
def checkPep(what,file,expected):
  success=True
  print "(*): checking posterior error probabilities for "+what+"..."
  processFile = open(file)
  psm_line=[]
  line = processFile.readline()
  # extract relevant lines from file
  while (line[0:5]!="PSMId"): # reading line by line, until "PSMId"
    line = processFile.readline()
  for i in range (1,101):
    if i==1: psm_line.append(processFile.readline()) # 1st psm
    elif i == 50: psm_line.append(processFile.readline()) # 50th psm
    elif i == 100: psm_line.append(processFile.readline()) # 100th psm
    else: processFile.readline() # skip line
  # extract peps from relevant lines
  pep=["","",""]
  i=0
  for i in range (0, len(pep)):
    tabsNumber = 0
    for j in range(0, len(psm_line[i])): # reading char by char, looking for 4th column
      if psm_line[i][j] == "\t":
        tabsNumber = tabsNumber + 1
      if tabsNumber == 3:
        pep[i] += psm_line[i][j] # pep of ith psm
  # check whether pep within 5% expected value
  for i in range (0,3):
    if float(pep[i])<expected[i]-(5*expected[i]/100) or float(pep[i])>expected[i]+(5*expected[i]/100):
      success=False
      print("...TEST FAILED: posterior error prob for "+what+" are outside the expected range")
      print("check "+file+" for details")
  return success

# performance increase when description of correct features option is enabled
def performanceD4On():
  success = True
  print "(*): checking performance with description of correct features option..."
  processFile = os.popen("grep \"New pi_0\" " + "/tmp/PERCOLATOR_D4on.txt")
  output = processFile.read()
  extracted_D4on = int(output[39:40])
  processFile = os.popen("grep \"New pi_0\" " + "/tmp/PERCOLATOR_psms.txt")
  output = processFile.read()
  extracted_D4off = int(output[39:40])
  if extracted_D4on < extracted_D4off:
    print("...TEST FAILED: percolator with -D 4 option performed worse than without it")
    print("check /tmp/PERCOLATOR_D4on.txt and /tmp/PERCOLATOR_D4off.txt for details" )
    success = False
  return success

psmFile=tmpDir+"/PERCOLATOR_psms.txt"
peptideFile=tmpDir+"/PERCOLATOR_peptides.txt"
proteinFile=tmpDir+"/PERCOLATOR_proteins.txt"
# number of significant psms within boundaries
success=checkNumberOfSignificant("psms",psmFile,283)
# number of significant peptrides within boundaries
success=checkNumberOfSignificant("peptides",peptideFile,221)
# number of significant proteins within boundaries
#success=checkNumberOfSignificant("proteins",proteinFile,153)
# psm: pi0 within boundaries
success=checkPi0("psms",psmFile,0.8912)
# peptides: pi0 within boundaries
success=checkPi0("peptides",peptideFile,0.9165)
# psm: pep within boundaries
expected=[2.61748e-13,3.26564e-09,7.28959e-08]
#success = checkPep("psms",psmFile, expected);
# peptide : pep within boundaries
expected=[4.47324e-14,3.52218e-09,1.7545e-07]
#success = checkPep("peptides",peptideFile, expected);
# performance increase with -D 4 option
#success = performanceD4On()

# if no errors were encountered, succeed
if success==True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)
