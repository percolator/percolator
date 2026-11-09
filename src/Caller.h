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
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Timer.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Enzyme.h"
#include "Globals.h"
#include "Normalizer.h"
#include "ProteinProbEstimator.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "SetHandler.h"
#include "TmpDir.h"

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
  Caller();
  virtual ~Caller();

  static string greeter();
  string extendedGreeter();
  bool parseOptions(int argc, char** argv);
  int run();

 protected:
  std::unique_ptr<Normalizer> pNorm_;
  int normalizerType_ = Normalizer::STDV;
  SanityCheck check_;
  std::unique_ptr<ProteinProbEstimator> protEstimator_;
  std::unique_ptr<Enzyme> enzyme_;
  // holds the concatenation of multiple input files, removed at the end
  std::unique_ptr<TmpDir> tmpDir_;

  // file input parameters
  bool tabInput_;
  bool readStdIn_;
  std::string inputFN_;
  std::vector<std::string> inputFNs_;
  std::string protEstimatorDecoyPrefix_;

  // file output parameters
  std::string tabOutputFN_;
  std::string weightOutputFN_;
  std::string psmResultFN_, peptideResultFN_, proteinResultFN_;
  std::string decoyPsmResultFN_, decoyPeptideResultFN_, decoyProteinResultFN_;
  bool outputRT_ = false;

  // report level parameters
  bool reportUniquePeptides_;
  bool targetDecoyCompetition_;
  bool useMixMax_;
  std::string inputSearchType_;

  // SVM / cross validation parameters
  double selectionFdr_, initialSelectionFdr_, testFdr_;
  unsigned int numIterations_, maxPSMs_, nestedXvalBins_, numThreads_;
  double selectedCpos_, selectedCneg_;
  bool classicTraining_ = false;
  bool reportEachIteration_, quickValidation_, trainBestPositive_,
      skipNormalizeScores_, analytics_, useResetAlgorithm_,
      useCompositionMatch_, useIrlsPep_, useInterpolatingPep_, usePavaPep_,
      useTdcPep_;

  // reporting parameters
  std::string call_;

  Timer timer;

  std::istream& getDataInStream(std::ifstream& fileStream);
  bool loadAndNormalizeData(std::istream& dataStream,
                            SetHandler& setHandler,
                            Scores& allScores);
  void calcAndOutputResult(Scores& allScores);

  void calculatePSMProb(Scores& allScores, bool uniquePeptideRun);
  void writeResults(Scores& allScores, bool unique, bool writeOutput);
  void calculateProteinProbabilities(Scores& allScores);
  void checkIsWritable(const std::string& filePath);
};

#endif /*CALLER_H_*/
