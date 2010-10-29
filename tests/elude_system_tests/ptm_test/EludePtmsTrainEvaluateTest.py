#!/usr/bin/python
# @Created by L. Moruz
# October 26th, 2010 
# Script that tests training and testing a model in Elude when ptms are present 

import os
import sys
sys.path.append('..')
import EludeUtilities as utility

# context format, test does not include rt 
# the alphabet of the train is same as the alphabet of the test
# remove in source and non enzymatic, check in source 
def TrainTestSameAlphabet():
  print "Running EludePtmsTrainEvaluateTest::TrainTestSameAlphabet..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_2.txt")
  out_file = os.path.join(path, "tmp.out")
  log_file = os.path.join(path, "tmp.log")
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file + " -o " + out_file
            + " -f -x -y -v 2 2> " + log_file)
  
  # check existence of output files 
  utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet", 
                              [out_file, log_file])

  # check content of the output file 
  lines = utility.loadFile(out_file)
  sorted(lines)
  if (len(lines) != 1191): 
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    number of peptides in the output file")
    return 
  
  peps = [lines[9].split("\t")[0], lines[1190].split("\t")[0]]
  rts = [float(lines[9].split("\t")[1]), float(lines[1190].split("\t")[1])]
  
  
  if (peps[0] != "K.KAEIGIS[unimod:21]MGSGTAVAK.S" or
      peps[1] != "R.TGRPALPVLYDNEEPFHIGQAK.V"):
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    peptides in the output file")
    return 
  
  if (abs(rts[0] - 19.39) > 5.0 or abs(rts[1] - 33.0) > 5.0):
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestSameAlphabet, incorrect \
    rts in the output file")
    return 
  
  # clean-up 
  utility.cleanUp([out_file, log_file])
  
  print "...TEST SUCCEEDED"

# context format, test includes rt 
# the alphabet of the train includes the alphabet of the test
# save and check the index 
def TrainTestDifferentAlphabet():
  print "Running EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_3.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  out_file = os.path.join(path, "tmp.out")
  log_file = os.path.join(path, "tmp.log")
  index_file = os.path.join(path, "tmp.index")
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file + " -o " + out_file 
            + " -r " + index_file + " -y -f -g -x -v 5 2> " + log_file)
  
  # check existence of output files 
  utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet", 
                              [out_file, log_file, index_file])
  
  # check if the index includes all the symbols 
  utility.checkIndex("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet", 
      index_file, [2, 10, 18, 20, 24], ["A", "K", "S[unimod:21]", "T[unimod:21]",
      "Y[unimod:21]"])
  
  # check content of the output file 
  utility.checkOutputFile("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet",
      out_file, 1410, [860, 1390, 1409], ["K.TVDAPILAAIK.K", 
      "K.HYGAGWLSMANAGADTNGSQFFITTVK.T", "R.IPQLQDVS[unimod:21]DFLKLLLLLLLLLLK.D"],
      [62.16, 118.92, 236.57], [58.51, 84.53, 83.68])

  # check performance 
  utility.checkPerformance("EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet",
                           log_file, (0.858, 0.877, 60.279))
  # clean-up 
  utility.cleanUp([out_file, log_file, index_file])
  print "...TEST SUCCEEDED"

# the alphabet of the train does not include the alphabet of the test
# context format, test includes rt 
def TrainTestInconsistentAlphabet():
  print "Running EludePtmsTrainEvaluateTest::TrainTestDifferentAlphabet..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(path, "tmp.log")
  out_file = os.path.join(path, "tmp.out")
  
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file  + 
            " -f -g -y -v 5 2> " + log_file)
  
  # check existence of output files
  utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet", 
                              [log_file])
  
  # check that no output file was generated 
  if os.path.isfile(out_file) :
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet, \
    Output file generated in spite of inconsistent alphabet.")
    return
  
  # clean-up 
  utility.cleanUp([log_file])
  print "...TEST SUCCEEDED" 
  
# context format, test includes rt, ignore ptms on
# the alphabet of the train does not include the alphabet of the test
# check performance measures 
def TrainTestInconsistentAlphabet():
  print "Running EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(path, "tmp.log")
  out_file = os.path.join(path, "tmp.out")
  
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file  + 
            " -f -g -y -v 5 2> " + log_file)
  
  # check existence of output files
  utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet", 
                              [log_file])
  
  # check that no output file was generated 
  if os.path.isfile(out_file) :
    utility.failTest("EludePtmsTrainEvaluateTest::TrainTestInconsistentAlphabet, \
    Output file generated in spite of inconsistent alphabet.")
    return
 
  # clean-up 
  utility.cleanUp([log_file])
  print "...TEST SUCCEEDED" 

# context format, test includes rt, ignore ptms on
# the alphabet of the train does not include the alphabet of the test
# check performance measures 
def TrainTestIgnorePtms():
  print "Running EludePtmsTrainEvaluateTest::TrainTestIgnorePtms..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_2.txt")
  test_file = os.path.join(data_folder, "test_3.txt")
  log_file = os.path.join(path, "tmp.log")
  out_file = os.path.join(path, "tmp.out")
  idx_file = os.path.join(path, "tmp.index")
  
  # run Elude (the flag -p indicates that ptms should be ignored)
  os.system(elude_path + " -t " + train_file + " -e " + test_file  + " -r " + idx_file 
            + " -o " + out_file  + " -f -g -p -v 5 2> " + log_file)
 
  # check existence of output files
  utility.checkFilesExistence("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms", 
                              [log_file, out_file])
  
  # check the index 
  # check the performance measure 
  utility.checkPerformance("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms",
                           log_file, (0.839, 0.854, 80.24))  
  
  utility.checkIndex("EludePtmsTrainEvaluateTest::TrainTestIgnorePtms",
      idx_file, [6, 19, 24], ["E[unimod:25]", "S[unimod:21]", "Y[unimod:21]"])
 
  # clean-up 
  utility.cleanUp([log_file, out_file, idx_file])
  print "...TEST SUCCEEDED" 
  
def main():
  TrainTestSameAlphabet()
  TrainTestDifferentAlphabet()
  TrainTestInconsistentAlphabet()
  TrainTestIgnorePtms()
  print ""  
  
if __name__ == '__main__':
  main()
 