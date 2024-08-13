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
  if "proteins" in what:
    output = getLine("Number of protein",file)
    extracted = output.split(": ")[1]
  else:
    output = getLine("New pi_0 estimate",file)
    extracted = output.split("yields ")[1].split(" target")[0]
  
  try:
    extracted = float(extracted)
  except ValueError:
    print("...TEST FAILED: could not read from " + file)
    extracted = 0
    success = False
  # check if number of significant within 10% of expected value 
  if abs(extracted - expected) > 0.1*expected: 
    print("...TEST FAILED: number of significant "+what+"=" + str(int(extracted)) + " is outside of desired range (" + str(expected) + " +/- 10%)")
    print("check "+file+" for details") 
    success = False
  return success

# check value of pi_0 is within 5% expected value
def checkPi0(what,file,expected):
  success = True
  print("(*): checking pi_0 estimate for "+what+"...")
  output = getLine("Selecting pi_0",file)
  extracted = output.split("pi_0=")[1]
  try:
    extracted = float(extracted)
  except ValueError:
    print("...TEST FAILED: could not read from " + file)
    extracted = 0
    success = False
    
  if abs(extracted - expected) > 0.05*expected:
    print("...TEST FAILED: "+what+" pi_0=" + str(extracted) + " is outside of desired range (" + str(expected) + " +/- 5%)")
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

if xmlSupport:
  print("- PERCOLATOR PIN XML FORMAT")

  psmFile = os.path.join(pathToOutputData,"PERCOLATOR_psms.txt")
  peptideFile = os.path.join(pathToOutputData,"PERCOLATOR_peptides.txt")
  proteinFile = os.path.join(pathToOutputData,"PERCOLATOR_proteins.txt")
  # proteinFileFido = os.path.join(pathToOutputData,"PERCOLATOR_proteins-fido.txt")

  # number of significant psms within boundaries
  success=checkNumberOfSignificant("psms",psmFile,1137) and success
  # number of significant peptrides within boundaries
  success=checkNumberOfSignificant("peptides",peptideFile,924) and success
  # number of significant proteins within boundaries for picked-protein
  success=checkNumberOfSignificant("proteins",proteinFile,340) and success
  # number of significant proteins within boundaries for fido (poorly calibrated)
  # success=checkNumberOfSignificant("proteins-fido",proteinFileFido,526) and success
  # psm: pi0 within boundaries
  success=checkPi0("psms",psmFile,0.8435) and success
  # peptides: pi0 within boundaries
  success=checkPi0("peptides",peptideFile,0.8655) and success
  # psm: pep within boundaries (old dataset)
  #expected=[2.61748e-13,3.26564e-09,7.28959e-08]
  #success = checkPep("psms",psmFile, expected);
  # peptide : pep within boundaries (old dataset)
  #expected=[4.47324e-14,3.52218e-09,1.7545e-07]
  #success = checkPep("peptides",peptideFile, expected);

print("- PERCOLATOR TAB FORMAT")

psmFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_psms.txt")
peptideFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_peptides.txt")
proteinFile = os.path.join(pathToOutputData,"PERCOLATOR_tab_proteins.txt")
# proteinFileFido = os.path.join(pathToOutputData,"PERCOLATOR_tab_proteins-fido.txt")

# number of significant psms within boundaries
success=checkNumberOfSignificant("psms",psmFile,1137) and success
# number of significant peptrides within boundaries
success=checkNumberOfSignificant("peptides",peptideFile,924) and success
# number of significant proteins within boundaries for picked-protein
success=checkNumberOfSignificant("proteins",proteinFile,340) and success
# number of significant proteins within boundaries for fido (poorly calibrated)
# success=checkNumberOfSignificant("proteins-fido",proteinFileFido,526) and success
# psm: pi0 within boundaries
success=checkPi0("psms",psmFile,0.8435) and success
# peptides: pi0 within boundaries
success=checkPi0("peptides",peptideFile,0.8655) and success

# if no errors were encountered, succeed
if success==True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)
