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

# context format, test includes rt 
# the alphabet of the train includes the alphabet of the test
# save and check the index 

# no context format, test includes rt 
# the alphabet of the train does not include the alphabet of the test

# no context format, test includes rt 
# the alphabet of the train does not include the alphabet of the test
# ignore ptms  on 

'''
# the training and test data are in format A.X.Y, test data does not include rt,
# the predictions and the model are saved
# check predictions and alphabet 
def TrainTestContextData():
  print "Running EludeTrainEvaluateTest::TrainTestContextData..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train.txt")
  test_file = os.path.join(data_folder, "test.txt")
  out_file = os.path.join(path, "tmp.out")
  model_file = os.path.join(path, "tmp.model")
  log_file = os.path.join(path, "tmp.log")
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file + " -o " + out_file
            + " -f -v 5 2> " + log_file)
  
  # check existence of output files 
  utility.checkFilesExistence("EludeTrainEvaluateTest::TrainTestContextData", 
                              [out_file, log_file])

  # check content of the output file 
  lines = utility.loadFile(out_file)
  sorted(lines)
  if len(lines) != 1255: 
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect number of \
    peptides in the output file")
    return   
  
  if lines[1093].split("\t")[0] != "K.IEGDVVVASAYSHELPR.Y":
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect peptide \
    in the output file")
    return 
  
  rt = float(lines[1093].split("\t")[1])
  if rt > 32 or rt < 22:
    utility.failTest("EludeTrainEvaluateTest::TrainTestContextData, incorrect predicted \
    rt in the output file")
    return   

  
  # clean-up 
  utility.cleanUp([out_file, log_file])
  
  print "...TEST SUCCEEDED"

# train and test when only peptide sequences are given and rts 
# I will not save the model, but I will check the performance figures from the log file
# also, check non enzymatic and in source 
def TrainTestNoContextData():
  print "Running EludeTrainEvaluateTest::TrainTestNoContextData..."
  path = os.path.dirname(sys.argv[0])
  
  elude_path = os.path.join(path, "./../../", "elude")
  data_folder = os.path.join(path, "./../../", "data/elude_test/standalone/")  
  train_file = os.path.join(data_folder, "train_1.txt")
  test_file = os.path.join(data_folder, "test_1.txt")
  log_file = os.path.join(path, "tmp.log")
  out_file = os.path.join(path, "tmp.out")
  insource_file = os.path.join(path, "in_source.txt")
  
  # run Elude  
  os.system(elude_path + " -t " + train_file + " -e " + test_file + " -o " 
            + out_file + " -i " + insource_file + " -g -x -u -y -v 5 2> " + log_file)
  
  # check existence of output files 
  utility.checkFilesExistence("EludeTrainEvaluateTest::TrainTestNoContextData", 
                              [insource_file, out_file, log_file])

  # check in source fragmentation
  utility.testFileContent(insource_file, 6, [3,4,5], 
                          ["YGASAGNVGDEGGVAPNIQTAEEALDLIVDAIK	105.947	test\n",
                           "IQFPHVADLLTSIQPPLTL	75.6913	train\n", 
                           "YQGYAEDVR	20.4553	test\n"])
  
  # check the number of lines in the output file   
  utility.testFileContent(out_file, 56, [], [])
  
  # check performance figures 
  utility.checkPerformance("EludeTrainEvaluateTest::TrainTestNoContextData", 
                           log_file, (0.87, 0.77, 73))
 
  # clean-up 
  utility.cleanUp([out_file, log_file, insource_file])
  
  print "...TEST SUCCEEDED"
'''

def main():
  #TrainTestContextData()
  print ""  
  
if __name__ == '__main__':
  main()
 