# Mattia Tomasoni - Percolator Project
# Script that tests the performances of percolator
# Parameters: none

import os
import sys
import re

pathToBinaries = "@pathToBinaries@"
pathToData = "@pathToData@"
pathToOutputData = "@pathToOutputData@"
xmlSupport = @xmlSupport@
success = True

print("PERCOLATOR PERFORMANCE")

def getLine(search_term,file):
  for line in open(file, 'r'):
    if line == None:
      print('no matches found')
      return ""
    if re.search(search_term, line):
      return line
    

# check number of psms/peptides with q-value < 0.01 is withing 5% from expected
def checkNumberOfSignificant(what,file,expected):
  success = True
  print("(*): checking number of significant "+what+" found...")
  if what=="proteins":
    output = getLine("The number of proteins idenfified at q-value = 0.01 is :",file)
  else:
    output = getLine("New pi_0 estimate",file)
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
  output = getLine("Selecting pi_0",file)
  extracted = float(output[15:20])
  if extracted<expected-(5*expected/100) or extracted>expected+(5*expected/100):
    print("...TEST FAILED: "+what+" pi_0=" + str(extracted) + " is outside of desired range")
    print("check "+file+" for details") 
    success = False
  return success

# check pep within 5% expected value
def checkPep(what,file,expected):
  success=True
  print("(*): checking posterior error probabilities for "+what+"...")
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
def performanceD4On(docFile, psmFile):
  success = True
  print("(*): checking performance with description of correct features option...")
  output = getLine("New pi_0", docFile)
  extracted_D4on = int(output[39:40])
  output = getLine("New pi_0", psmFile)
  extracted_D4off = int(output[39:40])
  if extracted_D4on < extracted_D4off:
    print("...TEST FAILED: percolator with -D 4 option performed worse than without it")
    print("check " + docFile + " and " + psmFile + " for details")
    success = False
  return success

if xmlSupport:
  print("- PERCOLATOR PIN XML FORMAT")

  psmFile = os.path.join(pathToOutputData,"PERCOLATOR_psms.txt")
  peptideFile = os.path.join(pathToOutputData,"PERCOLATOR_peptides.txt")
  proteinFile = os.path.join(pathToOutputData,"PERCOLATOR_proteins.txt")
  docFile = os.path.join(pathToOutputData,"PERCOLATOR_D4on.txt")

  # number of significant psms within boundaries
  success=checkNumberOfSignificant("psms",psmFile,292) and success
  # number of significant peptrides within boundaries
  success=checkNumberOfSignificant("peptides",peptideFile,221) and success
  # number of significant proteins within boundaries (old dataset)
  #success=checkNumberOfSignificant("proteins",proteinFile,153) and success
  # psm: pi0 within boundaries
  success=checkPi0("psms",psmFile,0.8912) and success
  # peptides: pi0 within boundaries
  success=checkPi0("peptides",peptideFile,0.9165) and success
  # psm: pep within boundaries (old dataset)
  #expected=[2.61748e-13,3.26564e-09,7.28959e-08]
  #success = checkPep("psms",psmFile, expected);
  # peptide : pep within boundaries (old dataset)
  #expected=[4.47324e-14,3.52218e-09,1.7545e-07]
  #success = checkPep("peptides",peptideFile, expected);
  # performance increase with -D 4 option
  success = performanceD4On(docFile, psmFile) and success

print("- PERCOLATOR TAB FORMAT")

psmFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_psms.txt")
peptideFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_peptides.txt")
proteinFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_proteins.txt")
docFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_D4on.txt")

# number of significant psms within boundaries (ubuntu=283, windows=301)
success=checkNumberOfSignificant("psms",psmFile,292) and success
# number of significant peptrides within boundaries (ubuntu=211, windows=225)
success=checkNumberOfSignificant("peptides",peptideFile,221) and success
# psm: pi0 within boundaries
success=checkPi0("psms",psmFile,0.8912) and success
# peptides: pi0 within boundaries
success=checkPi0("peptides",peptideFile,0.9165) and success
# performance increase with -D 4 option
success = performanceD4On(docFile, psmFile) and success

# if no errors were encountered, succeed
if success==True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)
