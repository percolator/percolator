#!/usr/bin/python
# @Created by L. Moruz
# October 25th, 2010 
# Script that tests loading a model in Elude when no ptms are present 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""

# trains a model, then loads the model and make sure that the predictions are the same 
# context data, test does not include rt
def SaveLoadContextData():
  print("Running EludeSaveLoadEvaluateTest::SaveLoadContextData...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train.txt")
  test_file = os.path.join(data_folder, "test.txt")
  out_file1 = os.path.join(out_path, "tmp1.out")
  out_file2 = os.path.join(out_path, "tmp2.out")
  model_file = os.path.join(out_path, "tmp.model")
  log_file = os.path.join(out_path, "tmp.log")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " + out_file1
            + " -s " + model_file + " -f -w -v 5 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeSaveLoadEvaluateTest::SaveLoadContextData", 
                              [out_file1, log_file, model_file]):
    utility.cleanUp([out_file1, log_file, model_file])
    exit(1)

  # check the alphabet in the model file 
  lines = utility.loadFile(model_file)
  if lines[len(lines) - 1] != "AA_alphabet 20 A C D E F G H I K L M N P Q R S T V W Y": 
    utility.failTest("EludeSaveLoadEvaluateTest::SaveLoadContextData, incorrect \
    alphabet in the model file")
    utility.cleanUp([out_file1, log_file, model_file])
    exit(1)
  
  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file + " -o " + out_file2
            + " -f -v 5 2> " + log_file)

  # check existence of the second output file
  if not utility.checkFilesExistence("EludeSaveLoadEvaluateTest::SaveLoadContextData", 
                              [out_file2]):
    utility.cleanUp([out_file1, log_file, model_file])
    exit(1)

  # check that the predictions match 
  lines1 = utility.loadFile(out_file1)
  lines2 = utility.loadFile(out_file2)
  sorted(lines1)
  sorted(lines2)
  if (len(lines1) != len(lines2)): 
    utility.failTest("EludeSaveLoadEvaluateTest::SaveLoadContextData, incorrect \
    number of peptides in the output file")
    utility.cleanUp([out_file1, out_file2, model_file, log_file])
    exit(1)

    
  fun = lambda x: (x.split("\t")[0], float(x.split("\t")[1]))
  for i in range(3, len(lines1), 100): 
    pep1, rt1 = fun(lines1[i])
    pep2, rt2 = fun(lines2[i])
    if (pep1 != pep2):
      utility.failTest("EludeSaveLoadEvaluateTest::SaveLoadContextData, incorrect \
      peptides in the output file")
      return 
    if (abs(rt1 - rt2) > 1.0): 
      utility.failTest("EludeSaveLoadEvaluateTest::SaveLoadContextData, incorrect \
      rts in the output file")
      return 

  # clean-up 
  utility.cleanUp([out_file1, out_file2, model_file, log_file])
 
  print("...TEST SUCCEEDED")
  
# trains a model, then loads the model and makes sure that the performance figures are the same 
# train includes rt, no context 
def SaveLoadNoContextData():
  print("Running EludeSaveLoadEvaluateTest::SaveLoadNoContextData...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_1.txt")
  test_file = os.path.join(data_folder, "test_1.txt")
  model_file = os.path.join(out_path, "tmp.model")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  
  # train and save the model
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file
            + " -s " + model_file + " -w -g -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeSaveLoadEvaluateTest::SaveLoadNoContextData", 
                              [log_file1, model_file]):
    utility.cleanUp([log_file1, model_file])
    exit(1)

  # run elude but load the model this time 
  os.system(pathToBinaries + " -l " + model_file + " -e " + test_file 
            + " -w -g -v 5 2> " + log_file2)

  if not utility.checkFilesExistence("EludeSaveLoadEvaluateTest::SaveLoadNoContextData", 
                              [log_file2]):
    utility.cleanUp([log_file1, model_file])
    exit(1)

  # check that the performance measures match 
  p1, s1, delta1 = utility.checkPerformance("", log_file1)
  p2, s2, delta2 = utility.checkPerformance("", log_file2)
  if abs(p1 - p2) > 0.1 or abs(s1 - s2) > 0.1 or abs(delta1 - delta2) > (delta1 * 0.05): 
    utility.failTest("EludeSaveLoadEvaluateTest::SaveLoadNoContextData, incorrect \
    performance figures")
    utility.cleanUp([model_file, log_file1, log_file2])
    exit(1)
  
  # clean-up 
  utility.cleanUp([model_file, log_file1, log_file2])
 
  print("...TEST SUCCEEDED")
  
def main():
  SaveLoadContextData()
  SaveLoadNoContextData()
  print("")

  
if __name__ == '__main__':
  main()
 
