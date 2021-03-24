#!/usr/bin/python
# @Created by L. Moruz
# October 15th, 2010 
# Script that tests training testing a model in Elude when no ptms are present 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""

# context data, test does not include rt
# the predictions are checked
def TrainTestContextData():
  print("Running EludeTrainEvaluateTest::TrainTestContextData...")
  print("!!!!!", pathToData, "\n", out_path, "\n", pathToBinaries, "\n!!!!!!!")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
 
  train_file = os.path.join(data_folder, "train.txt")
  test_file = os.path.join(data_folder, "test.txt")
  out_file = os.path.join(out_path, "tmp.out")
  log_file = os.path.join(out_path, "tmp.log")
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " + out_file
            + " -f -v 5 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeTrainEvaluateTest::TrainTestContextData", 
                              [out_file, log_file]):
    utility.cleanUp([out_file, log_file])
    exit(1)

  # check content of the output file 
  lines = utility.loadFile(out_file)
  sorted(lines)
  if len(lines) != 1255: 
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect number of \
    peptides in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1)
  
  if lines[1093].split("\t")[0] != "K.IEGDVVVASAYSHELPR.Y":
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect peptide \
    in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1) 
  
  rt = float(lines[1093].split("\t")[1])
  if rt > 32 or rt < 22:
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect predicted \
    rt in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1)
  
  # clean-up 
  utility.cleanUp([out_file, log_file])
  
  print("...TEST SUCCEEDED")

# train and test when only peptide sequences are given and rts 
# I will not save the model, but I will check the performance figures from the log file
# also, check non enzymatic and in source 
def TrainTestNoContextData():
  print("Running EludeTrainEvaluateTest::TrainTestNoContextData...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_1.txt")
  test_file = os.path.join(data_folder, "test_1.txt")
  log_file = os.path.join(out_path, "tmp.log")
  out_file = os.path.join(out_path, "tmp.out")
  insource_file = os.path.join(out_path, "in_source.txt")
  
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " 
            + out_file + " -i " + insource_file + " -g -x -u -y -v 5 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeTrainEvaluateTest::TrainTestNoContextData", 
                              [insource_file, out_file, log_file]):
    #utility.cleanUp([insource_file, out_file, log_file])
    exit(1)
  
  # check in source fragmentation
  if not utility.testFileContent(insource_file, 5, [3,4], 
                          ["IQFPHVADLLTSIQPPLTL	75.6913	train\n", 
                           "YQGYAEDVR	20.4553	test\n"]):
    utility.failTest("Incorrect content of the in-source fragmentation file")
    utility.cleanUp([insource_file, out_file, log_file])
    exit(1)
  
  # check the number of lines in the output file   
  if not utility.testFileContent(out_file, 56, [], []):
    utility.failTest("Incorrect number of lines in the output file")
    utility.cleanUp([insource_file, out_file, log_file])
    exit(1)
    
  # check performance figures 
  if utility.checkPerformance("EludeTrainEvaluateTest::TrainTestNoContextData", 
                           log_file, (0.87, 0.77, 73)) == (None, None, None):
    utility.cleanUp([insource_file, out_file, log_file])
    exit(1)

  # clean-up 
  utility.cleanUp([out_file, log_file, insource_file])
  
  print("...TEST SUCCEEDED")
  
def main():
  TrainTestContextData()
  TrainTestNoContextData()
  print(""  )
  
if __name__ == '__main__':
  main()
 
