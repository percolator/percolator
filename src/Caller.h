 /*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef CALLER_H_
#define CALLER_H_
#include "Timer.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <memory>
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif
#include "Globals.h"
#include "MyException.h"
#include "Option.h"
#include "SetHandler.h"
#include "DataSet.h"
#include "Scores.h"
#include "SanityCheck.h"
#include "Normalizer.h"
#include "ProteinProbEstimator.h"
#include "FidoInterface.h"
#include "PickedProteinInterface.h"
#include "XMLInterface.h"
#include "CrossValidation.h"
#include "Enzyme.h"

#define  NO_BOOST_DATE_TIME_INLINE
#include <boost/asio.hpp>
#include <boost/functional/hash_fwd.hpp>

#include <algorithm>


/*
* Main class that starts and controls the calculations.
*
* In addition to the calculation Caller also handles command line
* initialization, parsing input parameters, handling output of the
* calculation results.
* 
*/
class Caller {

 public:
  enum SetHandlerType {
    NORMAL = 1, SHUFFLED = -1, SHUFFLED_TEST = 2, SHUFFLED_THRESHOLD = 3
  };
  
  Caller();
  virtual ~Caller();
  
  static string greeter();
  string extendedGreeter();
  bool parseOptions(int argc, char **argv);    
  int run();
  
 protected:    
  Normalizer* pNorm_;
  SanityCheck* pCheck_;
  ProteinProbEstimator* protEstimator_;
  Enzyme* enzyme_;
  
  // file input parameters
  bool tabInput_;
  bool readStdIn_;
  std::string inputFN_;
  bool xmlSchemaValidation_;
  
  // file output parameters
  std::string tabOutputFN_, xmlOutputFN_, pepXMLOutputFN_;
  std::string weightOutputFN_;
  std::string psmResultFN_, peptideResultFN_, proteinResultFN_;
  std::string decoyPsmResultFN_, decoyPeptideResultFN_, decoyProteinResultFN_;
  bool xmlPrintDecoys_, xmlPrintExpMass_;
  
  // report level parameters
  bool reportUniquePeptides_;
  bool reportPepXML_;
  bool targetDecoyCompetition_;
  bool useMixMax_;
  std::string inputSearchType_;
  
  // SVM / cross validation parameters
  double selectionFdr_, initialSelectionFdr_, testFdr_;
  unsigned int numIterations_, maxPSMs_, nestedXvalBins_, numThreads_;
  double selectedCpos_, selectedCneg_;
  bool reportEachIteration_, quickValidation_, trainBestPositive_,
    skipNormalizeScores_;
  
  // reporting parameters
  std::string call_;

  Timer timer;
  
  std::istream& getDataInStream(std::ifstream& fileStream);
  bool loadAndNormalizeData(std::istream &dataStream, XMLInterface& xmlInterface, SetHandler& setHandler, Scores& allScores);
  void calcAndOutputResult(Scores& allScores, XMLInterface& xmlInterface);
  
  void calculatePSMProb(Scores& allScores, bool uniquePeptideRun);
  void calculateProteinProbabilities(Scores& allScores);
  void checkIsWritable(const std::string& filePath);
  
#ifdef CRUX
  virtual void processPsmScores(Scores& allScores) {}
  virtual void processPeptideScores(Scores& allScores) {}
  virtual void processProteinScores(ProteinProbEstimator* protEstimator) {}
#endif

bool detect_tab(std::string file_name)
{
    // open C++ stream to file
    std::ifstream file(file_name.c_str());
    // file not opened, return false
    if(!file.is_open()) return false;
    // read a line from the file       
    std::string wtf;
    std::istream &in= std::getline(file, wtf);
    // unable to read the line, return false
    if(!in) return false;
    // try to find a comma, return true if comma is found within the string
    return std::find(wtf.begin(), wtf.end(), '\t')!= wtf.end();
}
    
};



#endif /*CALLER_H_*/
