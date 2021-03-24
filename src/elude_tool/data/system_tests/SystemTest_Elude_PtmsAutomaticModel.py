#!/usr/bin/python
# @Created by L. Moruz
# October 27th, 2010 
# Script that tests automatic model selection in Elude when ptms are present 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""

# test add a model to library 
def AddModelLibrary():
  print("Running EludePtmsAutomaticModelTest::AddModelLibrary...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  log_file = os.path.join(out_path, "tmp.log")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -b " + out_path +  " -d -f -v 5 2> " 
            + log_file)
  
  model_file = os.path.join(out_path, "train_2.model")
  if not utility.checkFilesExistence("EludePtmsAutomaticModelTest::AddModelLibrary",
                              [log_file, model_file]):
    try:
      os.remove(log_file)
      os.remove(model_file)
      exit(1)
    except:
      exit(1)
    
  # clean-up 
  utility.cleanUp([log_file, model_file]) 
  print("...TEST SUCCEEDED")

# add different models to the library, make sure the suitable ones are chosen, no calibration 
# check output and performance, the ignore ptms flag is off
def AutomaticSelection():
  print("Running EludePtmsAutomaticModelTest::AutomaticSelection...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_files = map(lambda f: os.path.join(data_folder, f), ["train_2.txt", "train_3.txt"])
  test_file =  os.path.join(data_folder, "test_3.txt")
  log_file1 = os.path.join(out_path, "tmp1.log")
  out_file1 = os.path.join(out_path, "tmp1.out")
  log_file2 = os.path.join(out_path, "tmp2.log")
  out_file2 = os.path.join(out_path, "tmp2.out")
  lib_path = os.path.join(pathToData, "elude/calibrate_data/test_lib")
  
  # train 
  for t in train_files:
    os.system(pathToBinaries + " -t " + t + " -b " + lib_path +  " -d -f -v 1")
    
  expected_files = map(lambda f: os.path.join(lib_path, 
      os.path.split(f)[1].replace(".txt", ".model")), train_files)
  
  if not utility.checkFilesExistence("EludePtmsAutomaticModelTest::AutomaticSelection",
                              expected_files):
    utility.cleanUp(expected_files)
    exit(1)
    
  # run automatic model selection 
  os.system(pathToBinaries  + " -t " + train_files[1] + " -e " + test_file + " -o " 
            + out_file1 + " -b " + lib_path  +  " -j -f -g -w -a -v 5 2> " + log_file1)

  # just load the best model 
  os.system(pathToBinaries  + " -l " + expected_files[1] + " -e " + test_file + " -o " 
            + out_file2 +  " -f -g -v 5 2> " + log_file2)

  # check the performances 
  (p1, s1, d1) = utility.checkPerformance("", log_file1, None)   
  (p2, s2, d2) = utility.checkPerformance("", log_file2, None)   
  if (abs(p1 - p2) > 0.01 or abs(s1 - s2) > 0.01 or abs(d1 - d2) > 0.1):
    utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelection, incorrect \
    performance figures")
    utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
    exit(1)

  # check that the output are matching 
  lines1 = sorted(utility.loadFile(out_file1))
  lines2 = sorted(utility.loadFile(out_file2))
  if (len(lines1) != len(lines2)):
    utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelection, incorrect" +
        "number of lines in the output")
    utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
    exit(1)
  fun = lambda x: (x.split("\t")[0], float(x.split("\t")[1]), float(x.split("\t")[2]))
  for i in range(3, len(lines1), 100): 
    pep1, rt1, ort1 = fun(lines1[i])
    pep2, rt2, ort2 = fun(lines2[i])
    if (pep1 != pep2):
      utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelection, incorrect \
      peptides in the output file")
      utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
      exit(1)
    if (abs(rt1 - rt2) > 1.0): 
      utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelection, incorrect \
      predictions in the output file")      
      utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
      exit(1) 
    if (abs(ort1 - ort2) > 0.1): 
      utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelection, incorrect \
      observed rts in the output file")
      utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
      exit(1)
    
  # clean-up 
  utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
  
  print("...TEST SUCCEEDED")

# add different models to the library, make sure the suitable ones are chosen 
# check output and performance, the ignore ptms flag is on
def AutomaticSelectionIgnorePtms():
  print("Running EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_files = map(lambda f: os.path.join(data_folder, f), ["train_2.txt", "train_3.txt"])
  test_file =  os.path.join(data_folder, "test_3.txt")
  log_file1 = os.path.join(out_path, "tmp1.log")
  out_file1 = os.path.join(out_path, "tmp1.out")
  log_file2 = os.path.join(out_path, "tmp2.log")
  out_file2 = os.path.join(out_path, "tmp2.out")
  lib_path = os.path.join(pathToData, "elude/calibrate_data/test_lib")
  
  # train 
  for t in train_files:
    os.system(pathToBinaries + " -t " + t + " -b " + lib_path +  " -d -f -v 1")
    
  expected_files = map(lambda f: os.path.join(lib_path, 
      os.path.split(f)[1].replace(".txt", ".model")), train_files)
  
  if not utility.checkFilesExistence("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms",
                              expected_files):
    utility.cleanUp(expected_files)
    exit(1)
    
  # run automatic model selection 
  os.system(pathToBinaries  + " -t " + train_files[1] + " -e " + test_file + " -o " 
            + out_file1 + " -b " + lib_path  +  " -p -f -g -a -v 5 2> " + log_file1)

  # just load the best model 
  os.system(pathToBinaries  + " -l " + expected_files[1] + " -e " + test_file + " -o " 
            + out_file2 +  " -f -g -v 5 2> " + log_file2)

  if not utility.checkFilesExistence("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms",
                              [out_file1, log_file1, out_file2, log_file2]):
    utility.cleanUp(expected_files + [out_file1, log_file1, out_file2, log_file2])
    exit(1)
  
  # check the performances 
  (p1, s1, d1) = utility.checkPerformance("", log_file1, None)   
  (p2, s2, d2) = utility.checkPerformance("", log_file2, None)   
  if (abs(p1 - p2) > 0.01 or abs(s1 - s2) > 0.01):
    utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms, incorrect \
    performance figures")
    utility.cleanUp(expected_files + [out_file1, log_file1, out_file2, log_file2])
    exit(1)

  # check that the output are matching (only the peptides and the observed) 
  lines1 = sorted(utility.loadFile(out_file1))
  lines2 = sorted(utility.loadFile(out_file2))
  if (len(lines1) != len(lines2)):
    utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms, incorrect" +
        "number of lines in the output")
    utility.cleanUp(expected_files + [out_file1, log_file1, out_file2, log_file2])
    exit(1)

  fun = lambda x: (x.split("\t")[0], float(x.split("\t")[1]), float(x.split("\t")[2]))
  for i in range(3, len(lines1), 100): 
    pep1, rt1, ort1 = fun(lines1[i])
    pep2, rt2, ort2 = fun(lines2[i])
    if (pep1 != pep2):
      utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms, incorrect \
      peptides in the output file")
      utility.cleanUp(expected_files + [out_file1, log_file1, out_file2, log_file2])
      exit(1)
    if (abs(ort1 - ort2) > 0.1): 
      utility.failTest("EludePtmsAutomaticModelTest::AutomaticSelectionIgnorePtms, incorrect \
      observed rts in the output file")
      utility.cleanUp(expected_files + [out_file1, log_file1, out_file2, log_file2])
      exit(1)
  
  # clean-up 
  utility.cleanUp(expected_files + [log_file1, out_file1, log_file2, out_file2]) 
  print("...TEST SUCCEEDED")

def main():
  AddModelLibrary()
  AutomaticSelection()
  AutomaticSelectionIgnorePtms()
  print("")

if __name__ == '__main__':
  main()
