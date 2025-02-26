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
#include "CrossValidation.h"
#include "DataSet.h"
#include "Enzyme.h"
#include "Globals.h"
#include "MyException.h"
#include "Normalizer.h"
#include "Option.h"
#include "PickedProteinInterface.h"
#include "ProteinProbEstimator.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "SetHandler.h"
#include "ValidateTabFile.h"
#include "XMLInterface.h"

#define NO_BOOST_DATE_TIME_INLINE
#include <algorithm>
#include <boost/asio.hpp>
#include <boost/functional/hash_fwd.hpp>

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
        NORMAL = 1,
        SHUFFLED = -1,
        SHUFFLED_TEST = 2,
        SHUFFLED_THRESHOLD = 3
    };

    Caller();
    virtual ~Caller();

    static string greeter();
    string extendedGreeter();
    bool parseOptions(int argc, char** argv);
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
    std::vector<std::string> inputFNs_;
    bool xmlSchemaValidation_;
    std::string protEstimatorDecoyPrefix_;

    // file output parameters
    std::string tabOutputFN_, xmlOutputFN_, pepXMLOutputFN_;
    std::string weightOutputFN_;
    std::string psmResultFN_, peptideResultFN_, proteinResultFN_;
    std::string decoyPsmResultFN_, decoyPeptideResultFN_, decoyProteinResultFN_;
    bool xmlPrintDecoys_, xmlPrintExpMass_;
    bool outputRT_ = false;

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
        skipNormalizeScores_, analytics_, use_reset_alg_, use_composition_match_, use_spline_pep_, use_interpolating_pep_, use_pep_from_q_;

    // reporting parameters
    std::string call_;

    Timer timer;

    std::istream& getDataInStream(std::ifstream& fileStream);
    bool loadAndNormalizeData(std::istream& dataStream, XMLInterface& xmlInterface, SetHandler& setHandler, Scores& allScores);
    void calcAndOutputResult(Scores& allScores, XMLInterface& xmlInterface);

    void calculatePSMProb(Scores& allScores, bool uniquePeptideRun);
    void writeResults(Scores& allScores, bool unique, bool writeOutput);
    void calculateProteinProbabilities(Scores& allScores);
    void checkIsWritable(const std::string& filePath);
};

#endif /*CALLER_H_*/
