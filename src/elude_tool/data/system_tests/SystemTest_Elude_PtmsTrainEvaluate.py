#!/usr/bin/python
# @Created by L. Moruz
# October 26th, 2010 
# Script that tests training and testing a model in Elude when ptms are present 

import os
import sys
import SystemTest_Elude_Utilities as utility

pathToBinaries = os.path.join("", "elude")
pathToData = ""
out_path = ""

# context format, test does not include rt 
# the alphabet of the train is same as the alphabet of the test
# remove in source and non enzymatic, check in source 
def TrainTestSameAlphabet():
  print("Running EludePtmsTrainEvaluateTest::TrainTestSameAlphabet...")
   
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_2.txt")
  out_file = os.path.join(out_path, "tmp.out")
  log_file = os.path.join(out_path, "tmp.log")
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " + out_file
            + " -f -x -y -v 2 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet", 
                              [out_file, log_file]):
    utility.cleanUp([out_file, log_file])
    exit(1)

  # check content of the output file 
  lines = utility.loadFile(out_file)
  sorted(lines)
  if (len(lines) != 1191): 
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    number of peptides in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1)
  
  peps = [lines[9].split("\t")[0], lines[1190].split("\t")[0]]
  rts = [float(lines[9].split("\t")[1]), float(lines[1190].split("\t")[1])]  
  if (peps[0] != "K.KAEIGIS[unimod:21]MGSGTAVAK.S" or
      peps[1] != "R.TGRPALPVLYDNEEPFHIGQAK.V"):
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    peptides in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1)
    
  if (abs(rts[0] - 19.39) > 5.0 or abs(rts[1] - 33.0) > 5.0):
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    rts in the output file")
    utility.cleanUp([out_file, log_file])
    exit(1)
    
  # clean-up 
  utility.cleanUp([out_file, log_file])
  
  print("...TEST SUCCEEDED")

# context format, test includes rt which tests the in-source fragmentation functionality
# the alphabet of the train includes the alphabet of the test
# save and check the index 
def TrainTestDifferentAlphabet():
  print("Running EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_3.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  out_file = os.path.join(out_path, "tmp.out")
  log_file = os.path.join(out_path, "tmp.log")
  index_file = os.path.join(out_path, "tmp.index")
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file + " -o " + out_file 
            + " -r " + index_file + " -y -f -g -x -v 5 2> " + log_file)
  
  # check existence of output files 
  if not utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet", 
                              [out_file, log_file, index_file]):
    utility.cleanUp([out_file, log_file, index_file])
    exit(1)
  
  # check if the index includes all the symbols 
  if not utility.checkIndex("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet", 
      index_file, [2, 10, 18, 20, 24], ["A", "K", "S[unimod:21]", "T[unimod:21]",
      "Y[unimod:21]"]):
    utility.cleanUp([out_file, log_file, index_file])
    exit(1)
  
  # check content of the output file 
  if not utility.checkOutputFile("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet",
      out_file, 1410, [961, 1326, 1409], ["K.TVDAPILAAIK.K", 
      "K.HYGAGWLSMANAGADTNGSQFFITTVK.T", "R.IPQLQDVS[unimod:21]DFLKLLLLLLLLLLK.D"],
      [65.46, 105.63, 219.67], [58.51, 84.53, 83.68]):
    utility.cleanUp([out_file, log_file, index_file])
    exit(1)

  # check performance 
  if utility.checkPerformance("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet",
                           log_file, (0.858, 0.877, 60.279)) == (None, None, None):
    utility.cleanUp([out_file, log_file, index_file])
    exit(1)
    
  # clean-up 
  utility.cleanUp([out_file, log_file, index_file])
  print("...TEST SUCCEEDED")

# the alphabet of the train does not include the alphabet of the test
# context format, test includes rt 
def TrainTestInconsistentAlphabet():
  print("Running EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet...")

  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(out_path, "tmp.log")
  out_file = os.path.join(out_path, "tmp.out")
  
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file  + 
            " -f -g -y -v 5 2> " + log_file)
  
  # check existence of output files
  if not utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet", 
                              [log_file]):
    utility.cleanUp([log_file])
    exit(1)
  
  # check that no output file was generated 
  if os.path.isfile(out_file) :
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet, \
    Output file generated in spite of inconsistent alphabet.")
    utility.cleanUp([log_file, out_file])
    exit(1)
  
  # clean-up 
  utility.cleanUp([log_file])
  print("...TEST SUCCEEDED" )
  
# context format, test includes rt, ignore ptms on
# the alphabet of the train does not include the alphabet of the test
# check performance measures 
def TrainTestInconsistentAlphabet():
  print("Running EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet...")
  
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(out_path, "tmp.log")
  out_file = os.path.join(out_path, "tmp.out")
  
  # run Elude  
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file  + 
            " -f -g -y -v 5 2> " + log_file)
  
  # check existence of output files
  if not utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet", 
                              [log_file]):
    utility.cleanUp([log_file])
    exit(1)
    
  # check that no output file was generated 
  if os.path.isfile(out_file) :
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet, \
    Output file generated in spite of inconsistent alphabet.")
    utility.cleanUp([log_file, out_file])
    exit(1) 
    
  # clean-up 
  utility.cleanUp([log_file])
  print("...TEST SUCCEEDED" )

# context format, test includes rt, ignore ptms on
# the alphabet of the train does not include the alphabet of the test
# check performance measures 
def TrainTestIgnorePtms():
  print("Running EludePtmsTrainEvaluateTest::TrainTestIgnorePtms...")
 
  data_folder = os.path.join(pathToData, "elude/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(out_path, "tmp.log")
  out_file = os.path.join(out_path, "tmp.out")
  idx_file = os.path.join(out_path, "tmp.index")
  
  # run Elude (the flag -p indicates that ptms should be ignored)
  os.system(pathToBinaries + " -t " + train_file + " -e " + test_file  + " -r " + idx_file 
            + " -o " + out_file  + " -f -g -p -v 5 2> " + log_file)
 
  # check existence of output files
  if not utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms", 
                              [log_file, out_file, idx_file]):
    utility.cleanUp([log_file, out_file, idx_file])
    exit(1)
  
  # check the index 
  # check the performance measure 
  if utility.checkPerformance("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms",
                           log_file, (0.839, 0.854, 80.24)) == (None, None, None) or not utility.checkIndex("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms", idx_file, 
                             [6, 19, 24], ["E[unimod:25]", "S[unimod:21]", "Y[unimod:21]"]):
    utility.cleanUp([log_file, out_file, idx_file])
    exit(1)
 
  # clean-up 
  utility.cleanUp([log_file, out_file, idx_file])
  print("...TEST SUCCEEDED" )
  
def main():
  TrainTestSameAlphabet()
  TrainTestDifferentAlphabet()
  TrainTestInconsistentAlphabet()
  TrainTestIgnorePtms()
  print(""  )
  
if __name__ == '__main__':
  main()
 
