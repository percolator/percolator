#!/usr/bin/python
# @Created by L. Moruz
# October 25th, 2010 
# Script that tests automatic model selection in Elude 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""

# case 1: linear adjustment performed; the performances are compared to 
# when explicitly loading the best model 
# context data, rt in the test 
def AutomaticModelLinearAdjust():
  print("Running EludeAutomaticModelTest::AutomaticModelLinearAdjust...")
  #elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(pathToData, "elude/calibrate_data/")  
  lib_path = os.path.join(pathToData, "elude/calibrate_data/test_lib")  
  best_model = os.path.join(pathToData, "elude/calibrate_data/test_lib/Jupiter_60.model")  
  calibration_file = os.path.join(data_folder, "calibrate.txt")
  test_file = os.path.join(data_folder, "test.txt")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  
  # automatic model selection
  os.system(pathToBinaries + " -t " + calibration_file + " -e " + test_file + " -b " 
            + lib_path + " -w -g -f -a -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeAutomaticModelTest::AutomaticModelLinearAdjust", [log_file1]):
    exit(1)
  
  # run elude but load the best model this time 
  os.system(pathToBinaries + " -l " + best_model + " -e " + test_file + " -g -f " 
            + "-w -v 5 2> " + log_file2)

  # check existence of the second output file
  if not utility.checkFilesExistence("EludeAutomaticModelTest::AutomaticModelLinearAdjust", [log_file2]):
    utility.cleanUp([log_file1]) 
    exit(1)

  p1, s1, delta1 = utility.checkPerformance("", log_file1)
  p2, s2, delta2 = utility.checkPerformance("", log_file2)
  if abs(p1 - p2) > 0.1 or abs(s1 - s2) > 0.1: 
    utility.failTest("EludeAutomaticModelTest::AutomaticModelLinearAdjust, incorrect \
    performance figures")
    utility.cleanUp([log_file1, log_file2])
    exit(1)
  
  # clean-up 
  utility.cleanUp([log_file1, log_file2])

  print("...TEST SUCCEEDED")

# case 2: no linear adjustment performed; the performances are compared to 
# when explicitly loading the best model 
# context, test includes rt
def AutomaticModelNoLinearAdjust():
  print("Running EludeAutomaticModelTest::AutomaticModelNoLinearAdjust...")
  
  data_folder = os.path.join(pathToData, "elude/calibrate_data/")  
  lib_path = os.path.join(pathToData, "elude/calibrate_data/test_lib")  
  best_model = os.path.join(pathToData, "elude/calibrate_data/test_lib/yeast_20.model")  
  calibration_file = os.path.join(data_folder, "calibrate_1.txt")
  test_file = os.path.join(data_folder, "test_1.txt")
  log_file1 = os.path.join(out_path, "tmp1.log")
  log_file2 = os.path.join(out_path, "tmp2.log")
  
  # automatic model selection
  os.system(pathToBinaries + " -t " + calibration_file + " -e " + test_file + " -b " 
            + lib_path + " -w -x -y -j -g -f -a -v 5 2> " + log_file1)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludeAutomaticModelTest::AutomaticModelLinearNoAdjust", 
                              [log_file1]):
    exit(1)

  # run elude but load the best model this time 
  os.system(pathToBinaries + " -l " + best_model + " -e " + test_file + " -g -f " 
            + "-w -v 5 2> " + log_file2)

  # check existence of the second output file
  if not utility.checkFilesExistence("EludeAutomaticModelTest::AutomaticModelNoLinearAdjust", 
                              [log_file2]):
    utility.cleanUp([log_file1])
    exit(1)

  p1, s1, delta1 = utility.checkPerformance("", log_file1)
  p2, s2, delta2 = utility.checkPerformance("", log_file2)
  if abs(p1 - p2) > 0.1 or abs(s1 - s2) > 0.1 or abs(delta1 - delta2) > 1.0: 
    utility.failTest("EludeAutomaticModelTest::AutomaticModelNoLinearAdjust, incorrect \
    performance figures")
    utility.cleanUp([log_file1, log_file2])
    exit(1)
  
  # clean-up 
  utility.cleanUp([log_file1, log_file2])

  print("...TEST SUCCEEDED")

def main():
  AutomaticModelLinearAdjust()
  AutomaticModelNoLinearAdjust()
  print("")

if __name__ == '__main__':
  main()
  
