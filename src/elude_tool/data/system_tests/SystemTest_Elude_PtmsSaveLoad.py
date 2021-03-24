#!/usr/bin/python
# @Created by L. Moruz
# October 26th, 2010 
# Script that tests saving and loading a model in Elude when ptms are present 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""


# the model is trained on data with the same ptms as the test
def SaveLoadSamePtms():
  print("Running EludePtmsSaveLoadTest::SaveLoadSamePtms...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_2.txt")
  out_file1 = os.path.join(out_path, "tmp1.out")
  out_file2 = os.path.join(out_path, "tmp2.out")
  model_file = os.path.join(out_path, "tmp.model")
  log_file = os.path.join(out_path, "tmp.log")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " + out_file1
            + " -s " + model_file + " -f -v 5 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadSamePtms", 
                              [out_file1, log_file, model_file]):
    utility.cleanUp([out_file1, log_file, model_file])
    exit(1)

  # check the alphabet in the model file 
  lines = utility.loadFile(model_file)
  if lines[len(lines) - 1] != "AA_alphabet 23 A C D E E[unimod:25] F G H I K L M N P Q R S S[unimod:21] T V W Y Y[unimod:21]": 
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadSamePtms, incorrect \
    alphabet in the model file")
    utility.cleanUp([out_file1, log_file, model_file])
    exit(1)
  
  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file + " -o " + out_file2
            + " -f -v 5 2> " + log_file)
  
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadSamePtms", 
                              [out_file2]):
    utility.cleanUp([out_file1, log_file, model_file, out_file2])
    exit(1)
  
  # check that the predictions match  
  lines1 = utility.loadFile(out_file1)
  lines2 = utility.loadFile(out_file2)
  sorted(lines1)
  sorted(lines2)
  if (len(lines1) != len(lines2)): 
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadSamePtms, incorrect \
    number of peptides in the output file")
    utility.cleanUp([out_file1, out_file2, model_file, log_file]) 
    exit(1)
    
  fun = lambda x: (x.split("\t")[0], float(x.split("\t")[1]))
  for i in range(3, len(lines1), 100): 
    pep1, rt1 = fun(lines1[i])
    pep2, rt2 = fun(lines2[i])
    if (pep1 != pep2):
      utility.failTest("EludePtmsSaveLoadTest::SaveLoadSamePtms, incorrect \
      peptides in the output file")
      utility.cleanUp([out_file1, out_file2, model_file, log_file]) 
      exit(1)
    if (abs(rt1 - rt2) > 1.0): 
      utility.failTest("EludePtmsSaveLoadTest::SaveLoadSamePtms, incorrect \
      rts in the output file")
      utility.cleanUp([out_file1, out_file2, model_file, log_file]) 
      exit(1)

  # clean-up 
  utility.cleanUp([out_file1, out_file2, model_file, log_file]) 
  print("...TEST SUCCEEDED")
  
  
# the model is trained on data including the ptms in the test 
def SaveLoadDifferentPtms():
  print("Running EludePtmsSaveLoadTest::SaveLoadDifferentPtms...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_3.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  model_file = os.path.join(out_path, "tmp.model")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -s " 
            + model_file + " -w -f -g -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadDifferentPtms", 
                              [log_file1, model_file]):
    utility.cleanUp([log_file1, model_file])
    exit(1)

  # check the alphabet in the model file 
  lines = utility.loadFile(model_file)
  if lines[len(lines) - 1] != "AA_alphabet 23 A C D E F G H I K L M N P Q R S S[unimod:21] T T[unimod:21] V W Y Y[unimod:21]": 
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadDifferentPtms, incorrect \
    alphabet in the model file")
    utility.cleanUp([log_file1, model_file])
    exit(1)
  
  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file + " -w -f -g -v 5 2> " 
            + log_file2)

  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadDifferentPtms", 
                                     [log_file2]):
    utility.cleanUp([log_file1, model_file, log_file2])
    exit(1)  
  
  # check that the performances are the same 
  (p1, s1, d1) = utility.checkPerformance("", log_file1, None)   
  (p2, s2, d2) = utility.checkPerformance("", log_file2, None)   
  if (abs(p1 - p2) > 0.01 or abs(s1 - s2) > 0.01 or abs(d1 - d2) > 0.1):
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadDifferentPtms, incorrect \
    performance figures")
    utility.cleanUp([model_file, log_file1, log_file2]) 
    exit(1)
  
  # clean-up 
  utility.cleanUp([model_file, log_file1, log_file2]) 
  print("...TEST SUCCEEDED"  )
  
# the test includes ptms that were not present in the training set 
# the model is trained on data including the ptms in the test 
def SaveLoadInconsistentPtms():
  print("Running EludePtmsSaveLoadTest::SaveLoadInconsistentPtms...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  model_file = os.path.join(out_path, "tmp.model")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  out_file2 = os.path.join(out_path, "tmp2.out")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -s " 
            + model_file + " -f -g -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadInconsistentPtms", 
                              [log_file1, model_file]):
    utility.cleanUp([log_file1, model_file])
    exit(1)
  
  # check the alphabet in the model file 
  lines = utility.loadFile(model_file)
  if lines[len(lines) - 1] != "AA_alphabet 23 A C D E E[unimod:25] F G H I K L M N P Q R S S[unimod:21] T V W Y Y[unimod:21]": 
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadInconsistentPtms, incorrect \
    alphabet in the model file")
    utility.cleanUp([log_file1, model_file])
    exit(1)
  
  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file + " -o " + out_file2 
            + " -f -g -v 5 2> " + log_file2)
  
  if os.path.isfile(out_file2) :
    failTest("EludePtmsSaveLoadTest::SaveLoadInconsistentPtms, output file incorrectly" 
             + "generated")
    utility.cleanUp([log_file1, model_file, out_file2, log_file2])
    exit(1)
  
  # clean-up 
  utility.cleanUp([model_file, log_file1, log_file2]) 
  print("...TEST SUCCEEDED"  )

# the test includes ptms that were not present in the training set
# but the ignore ptms flag is set 
def SaveLoadIgnorePtms():
  print("Running EludePtmsSaveLoadTest::SaveLoadIgnorePtms...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  model_file = os.path.join(out_path, "tmp.model")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  out_file1 = os.path.join(out_path, "tmp1.out")
  out_file2 = os.path.join(out_path, "tmp2.out")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -s " + model_file 
            + " -o " + out_file1 + " -f -g -p -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadInconsistentPtms", 
                              [log_file1, model_file, out_file1]):
    utility.cleanUp([log_file1, model_file, out_file1])    
    exit(1)
   
  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file + " -o " + out_file2 
            + " -f -g -p -v 5 2> " + log_file2)
 
  if not utility.checkFilesExistence("EludePtmsSaveLoadTest::SaveLoadInconsistentPtms", 
                              [log_file2, out_file2]):
    utility.cleanUp([log_file1, model_file, out_file1, log_file2, out_file2])    
    exit(1)
  
  (p1, s1, d1) = utility.checkPerformance("", log_file1, None)   
  (p2, s2, d2) = utility.checkPerformance("", log_file2, None)   
  if (abs(p1 - p2) > 0.01 or abs(s1 - s2) > 0.01 or abs(d1 - d2) > 0.1):
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadIgnorePtms, incorrect \
    performance figures")
    utility.cleanUp([log_file1, model_file, out_file1, log_file2, out_file2])    
    exit(1)
  
   # check that the predictions match  
  lines1 = utility.loadFile(out_file1)
  lines2 = utility.loadFile(out_file2)
  sorted(lines1)
  sorted(lines2)
  if (len(lines1) != len(lines2)): 
    utility.failTest("EludePtmsSaveLoadTest::SaveLoadIgnorePtms, incorrect \
    number of peptides in the output file")
    utility.cleanUp([log_file1, model_file, out_file1, log_file2, out_file2])    
    exit(1)
    
  fun = lambda x: (x.split("\t")[0], float(x.split("\t")[1]))
  for i in range(3, len(lines1), 100): 
    pep1, rt1 = fun(lines1[i])
    pep2, rt2 = fun(lines2[i])
    if (pep1 != pep2):
      utility.failTest("EludePtmsSaveLoadTest::SaveLoadIgnorePtms, incorrect \
      peptides in the output file")
      utility.cleanUp([log_file1, model_file, out_file1, log_file2, out_file2])    
      exit(1)
    if (abs(rt1 - rt2) > 1.0): 
      utility.failTest("EludePtmsSaveLoadTest::SaveLoadIgnorePtms, incorrect \
      rts in the output file")
      utility.cleanUp([log_file1, model_file, out_file1, log_file2, out_file2])    
      exit(1)

  # clean-up 
  utility.cleanUp([out_file1, out_file2, model_file, log_file1, log_file2]) 
  print("...TEST SUCCEEDED"  )


def main():
  SaveLoadSamePtms()
  SaveLoadDifferentPtms()
  SaveLoadInconsistentPtms()
  SaveLoadIgnorePtms()
  
if __name__ == '__main__':
  main()
