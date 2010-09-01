/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *****************************************************************************/
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <ctime>
#include <sys/stat.h>
#include "Option.h"
#include "Globals.h"
#include "Elude.h"
#include "PSMDescription.h"
#include "LTSRegression.h"

using namespace std;

// path to the library
string RTPredictor::libPath =
    "models/";

// difference in hydrophobicity between parend and child when detecting CID fragments
float RTPredictor::diff_hydrophobicity = 10.0;

// how many peptides included in the time window reported
float RTPredictor::fractionPeptides = 0.95;

// used to sort a vector of pairs PSMDescription, bool
struct myPair {
    bool operator()(pair<pair<PSMDescription, bool> , bool> psm1, pair<
        pair<PSMDescription, bool> , bool> psm2) {
      return (psm1.first.first.getRetentionTime()
          < psm2.first.first.getRetentionTime());
    }
} mypair;

// check if a psm is CID fragment and in the training data
bool decayTrain(pair<pair<PSMDescription, bool> , bool> psm) {
  if (psm.second && (psm.first.second)) {
    return true;
  }
  return false;
}

// check if a psm is CID fragment
bool decay(pair<pair<PSMDescription, bool> , bool> psm) {
  if (psm.second) {
    return true;
  }
  return false;
}

// check if a psm is in the training data
bool inTrain(pair<pair<PSMDescription, bool> , bool> psm) {
  if (psm.first.second) {
    return true;
  }
  return false;
}

// get the psm
PSMDescription getPSM(pair<pair<PSMDescription, bool> , bool> psm) {
  return psm.first.first;
}

// compare psms according to retention time
bool comparePsmsRT(PSMDescription psm1, PSMDescription psm2) {
  return psm1.getRetentionTime() < psm2.getRetentionTime();
}

// compare psms according to predicted retention time
bool comparePsmsPRT(pair<PSMDescription, int> psm1, pair<PSMDescription,
    int> psm2) {
  double prt1, prt2;
  prt1 = psm1.first.getPredictedRetentionTime();
  prt2 = psm2.first.getPredictedRetentionTime();
  return prt1 < prt2;
}

// compare psms according to absolute value of difference between predicted and observed
// compare psms according to predicted retention time
bool comparePsmsDeltaRT(PSMDescription psm1, PSMDescription psm2) {
  double deltaRT1, deltaRT2;
  deltaRT1 = psm1.getPredictedRetentionTime() - psm1.getRetentionTime();
  deltaRT2 = psm2.getPredictedRetentionTime() - psm2.getRetentionTime();
  return deltaRT1 < deltaRT2;
}

// constructor
RTPredictor::RTPredictor() :
  trainFile(""), testFile(""), outputFile(""), saveModelFile(""),
      loadModelFile(""), trainRetentionFeatures(NULL),
      testRetentionFeatures(NULL), theNormalizer(NULL),
      decayingPeptidesFile(""), hydrophobicityFile(""),
      removeRedundant(false), autoModelSelection(false),
      linearAdjust(true), addLibModel(false), removeDecaying(false),
      testIncludesRT(true), removeNEnzymatic(false), theEnzyme(NULL) {
  RTModel::setDoKlammer(false);
  Normalizer::setType(Normalizer::UNI);
  Enzyme::setEnzyme(Enzyme::TRYPSIN);
}

// destructor
RTPredictor::~RTPredictor() {
  if (trainRetentionFeatures) {
    // delete the retention feature table
    delete[] trainRetentionFeatures;
    trainRetentionFeatures = NULL;
  }
  if (testRetentionFeatures) {
    // delete the retention feature table
    delete[] testRetentionFeatures;
    testRetentionFeatures = NULL;
  }
  if (theNormalizer) {
    delete theNormalizer;
  }
}

// introductory message
string RTPredictor::greeter() {
  ostringstream oss;
  oss << "Elude version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under MIT License" << endl;
  oss
      << "Written by Lukas Käll (lukas.kall@cbr.su.se) and Luminita Moruz (lumi@sbc.su.se)"
      << endl;
  oss << "Center for Biomembrane Reasearch." << endl;
  oss << "Dept. of Biochemistry, Stockholm University, Stockholm." << endl;
  oss << "Usage:" << endl;
  oss << "   elude [options]" << endl << endl;
  return oss.str();
}

// parse the command line
bool RTPredictor::parseOptions(int argc, char** argv) {
  string enzyme;
  ostringstream intro;
  intro << greeter() << endl;
  CommandLineParser cmd(intro.str());
  // define available options
  cmd.defineOption("v",
                   "verbose",
                   "Set verbosity of output: 0 = no processing info, 5 = all, default is 2",
                   "level");
  cmd.defineOption("n",
                   "no_grid",
                   "Specifies that no calibration of parameters should be carried out",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("f",
                   "fine_grid",
                   "Specifies that the calibration of parameters will be carried out using a coarse grid followed by a fine grid",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("p",
                   "plain_evaluation",
                   "The performance of a model is estimated by training the model on 3/4 of training data and testing on the 1/4 left.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("t",
                   "train",
                   "Specifies the file including the training data",
                   "filename");
  cmd.defineOption("c",
                   "calibration",
                   "File storing information about calibration steps",
                   "filename");
  cmd.defineOption("e",
                   "evaluate",
                   "Specifies the file including the test data",
                   "filename");
  cmd.defineOption("s",
                   "save-model",
                   "Specifies the file in which the model will be saved",
                   "filename");
  cmd.defineOption("m",
                   "load-model",
                   "Specifies a file including a SVR model to be loaded",
                   "filename");
  cmd.defineOption("o",
                   "out",
                   "Create an output file including information regarding the current run",
                   "filename");
  cmd.defineOption("u",
                   "unique",
                   "Remove all redundant peptides from the test set",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("k",
                   "k_fold",
                   "Specify the number of folds for cross validation",
                   "value");
  cmd.defineOption("r",
                   "rt_feat",
                   "Specify a number between 2^0 and 2^12-1 to select the rt feature groups used ",
                   "value");
  cmd.defineOption("a",
                   "auto",
                   "Specifies that the SVR model employed is the one that best matches the training data",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("d",
                   "append",
                   "Append current model to library ",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("j",
                   "no_linear_adjust",
                   "Specifies that the model will not be linearly adjusted ",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("y",
                   "no_in_source",
                   "Specifies that in source fragments should be removed from the test set",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("i",
                   "save-decay-peptides",
                   "Specifies the file in which the CID fragments will be stored",
                   "filename");
  cmd.defineOption("z",
                   "enzyme",
                   "Specifies the enzyme used for digestion. Possible values {NO_ENZYME,TRYPSIN,CHYMOTRYPSIN,ELASTASE}. By default TRYPSIN",
                   "value");
  cmd.defineOption("x",
                   "remove-non-enzymatic",
                   "All non enzymatic peptides should be removed from both train and test ",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("g",
                   "hydrophobicity",
                   "Specifies a file where the hydrophobicity index trained will be saved",
                   "filename");
  cmd.defineOption("b",
                   "lts_coverage",
                   "Specifies the fraction of data used in calibrating a model via LTS",
                   "value");
  // EXPERIMENTAL; specify the slack penalty that is to be used
  //cmd.defineOption("b", "slack_penalty", "Specify the slack penalty for the index SVR ", "value");
  // parse command line
  cmd.parseArgs(argc, argv);
  // process options
  if (cmd.optionSet("c")) {
    model.setCalibrationFile(cmd.options["c"]);
  }
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (cmd.optionSet("t")) {
    trainFile = cmd.options["t"];
  }
  if (cmd.optionSet("e")) {
    testFile = cmd.options["e"];
  }
  if (cmd.optionSet("s")) {
    saveModelFile = cmd.options["s"];
  }
  if (cmd.optionSet("m")) {
    loadModelFile = cmd.options["m"];
  }
  if (cmd.optionSet("o")) {
    outputFile = cmd.options["o"];
  }
  if (cmd.optionSet("n")) {
    model.setGridType(NO_GRID);
  }
  if (cmd.optionSet("f")) {
    model.setGridType(FINE_GRID);
  }
  if (cmd.optionSet("p")) {
    model.setEvaluationType(SIMPLE_EVAL);
  }
  if (cmd.optionSet("u")) {
    removeRedundant = true;
  }
  if (cmd.optionSet("k")) {
    model.setK(cmd.getInt("k", 1, 100));
  }
  if (cmd.optionSet("r")) {
    model.setSelectFeatures(cmd.getInt("r", 1, 1000));
  }
  if (cmd.optionSet("a")) {
    autoModelSelection = true;
  }
  if (cmd.optionSet("d")) {
    addLibModel = true;
  }
  if (cmd.optionSet("j")) {
    linearAdjust = false;
  }
  if (cmd.optionSet("y")) {
    removeDecaying = true;
  }
  if (cmd.optionSet("i")) {
    decayingPeptidesFile = cmd.options["i"];
  }
  if (cmd.optionSet("x")) {
    removeNEnzymatic = true;
  }
  if (cmd.optionSet("z")) {
    enzyme = cmd.options["z"];
  }
  if (!enzyme.empty()) {
    if ((enzyme.compare("CHYMOTRYPSIN") == 0)
        || (enzyme.compare("chymotrypsin") == 0)) {
      theEnzyme->setEnzyme(Enzyme::CHYMOTRYPSIN);
    } else if ((enzyme.compare("ELASTASE") == 0)
        || (enzyme.compare("elastase") == 0)) {
      theEnzyme->setEnzyme(Enzyme::ELASTASE);
    } else if ((enzyme.compare("TRYPSIN") == 0)
        || (enzyme.compare("trypsin") == 0)) {
      theEnzyme->setEnzyme(Enzyme::TRYPSIN);
    } else {
      theEnzyme->setEnzyme(Enzyme::NO_ENZYME);
    }
  }
  if (cmd.optionSet("g")) {
    hydrophobicityFile = cmd.options["g"];
  }
  if (cmd.optionSet("b")) {
    double coverage = cmd.getDouble("b", 0.0, 1.0);
    LTSRegression::setCoverage(coverage);
  }
  // EXPERIMENTAL; specify the slack penalty that is to be used when generating our hydrophobicity index
  //if (cmd.optionSet("b"))
  // C = cmd.getDouble("b", -100000.0, 100000.0);
  return true;
}

// print information about the current run
void RTPredictor::printRunInformation() {
  ostringstream oss;
  oss << "----------------------------------------" << endl;
  oss << "Running Elude v" << VERSION << " with the following options: "
      << endl;
  if (!trainFile.empty()) {
    oss << " *Training data: " << trainFile << endl;
  }
  if (!testFile.empty()) {
    oss << " *Test data: " << testFile << endl;
  }
  if (!loadModelFile.empty()) {
    oss << " *Loading model from file: " << loadModelFile << endl;
  }
  if (!saveModelFile.empty()) {
    oss << " *Saving model to file: " << saveModelFile << endl;
  }
  if (!hydrophobicityFile.empty()) {
    oss << " *Saving hydrophobicity model to file: " << hydrophobicityFile
        << endl;
  }
  if (autoModelSelection) {
    oss << " *The SVR model will be selected automatically" << endl;
    if (linearAdjust) {
      oss << " *Linear fitting will be carried out " << endl;
    } else {
      oss << " *No linear fitting will be carried out " << endl;
    }
  }
  if (!outputFile.empty()) {
    oss << " *Output file: " << outputFile << endl;
  }
  if (addLibModel) {
    oss << " *Current model will be added to the library" << endl;
  }
  if (removeRedundant) {
    oss << " *Redundant peptides will be removed from the test set"
        << endl;
  } else {
    oss << " *Redundant peptides are allowed in the test set" << endl;
  }
  if (removeDecaying) {
    oss << " *In source fragments will be removed from the test set"
        << endl;
  }
  if (!decayingPeptidesFile.empty()) {
    oss << " *In source fragments are saved to file "
        << decayingPeptidesFile << endl;
  }
  if (removeNEnzymatic) {
    oss
        << " *Non-enzymatic peptides are removed from both train and test set"
        << endl;
  }
  if (!trainFile.empty() && !autoModelSelection) {
    oss << " *Calibration of parameters: " << model.getGridType() << endl;
    oss << " *Evaluation type: " << model.getEvaluationType() << endl;
  }
  oss << " *Enzyme: " << theEnzyme->getStringEnzyme() << endl;
  if (!trainFile.empty()) {
    model.printFeaturesInUse(oss);
  }
  oss << "-----------------------------------------" << endl;
  if (VERB > 2) {
    cerr << oss.str();
  }
}

/*
 * LOAD AND PROCESS TRAINING DATA
 */
// load training data
void RTPredictor::loadTrainFile() {
  // retention time
  double retTime;
  string peptide;
  // opens the file for reading
  ifstream in(trainFile.c_str(), ios::in);
  if (VERB > 2) {
    cerr << endl << "Loading " << trainFile << "..." << endl;
  }
  while (in >> peptide >> retTime) {
    trainPsms.push_back(PSMDescription(peptide, retTime));
  }
  if (VERB > 2) {
    cerr << trainPsms.size() << " peptides were loaded." << endl;
  }
}

/*
 // process train data (load input file, remove redundant peptides, in source CID fragments and non tryptic peptides)
 void RTPredictor::processTrainData()
 {
 vector< pair<pair<PSMDescription, bool>,bool> >  psmPairs;

 // load input data
 loadTrainFile();

 // remove redundant peptides
 removeRedundantPeptides(trainPsms);

 if (!testFile.empty())
 {
 // load test data
 loadTestFile();

 // remove redundant peptides
 if (removeRedundant)
 removeRedundantPeptides(testPsms);
 else
 sort(testPsms.begin(), testPsms.end(), less<PSMDescription>());

 // remove from the train set the peptides featuring both in train and test
 if (VERB > 2 && !autoModelSelection)
 cerr << endl << "Removing from the train set the peptides featuring in both train and test sets..." << endl;

 if (!autoModelSelection)
 {
 vector<PSMDescription>::iterator it;
 int initialSize = trainPsms.size();
 for(it = testPsms.begin(); it != testPsms.end(); ++it)
 trainPsms.erase(remove(trainPsms.begin(), trainPsms.end(), (*it)), trainPsms.end());

 if (VERB > 2)
 {
 cerr << (initialSize - trainPsms.size())<< " peptides were removed." << endl;
 cerr << "Train set includes now " << trainPsms.size() << " peptides." << endl;
 }
 }

 // to detect CID fragments
 if (testIncludesRT)
 addToPairVector(testPsms, false, psmPairs);
 }

 addToPairVector(trainPsms, true, psmPairs);

 removeSourceCIDs(psmPairs);

 // divide in test and train
 if (testPsms.empty() || !testIncludesRT)
 {
 trainPsms.resize(psmPairs.size());
 transform(psmPairs.begin(), psmPairs.end(), trainPsms.begin(), getPSM);
 if (VERB > 2)
 cerr << "Train set includes now " << trainPsms.size() << " peptides." <<  endl;
 }
 else
 {
 int sizeTest = testPsms.size();
 vector< pair<pair<PSMDescription, bool>,bool> >::iterator it;
 it = partition(psmPairs.begin(), psmPairs.end(), inTrain);
 trainPsms.resize(distance(psmPairs.begin(), it));
 testPsms.resize(distance(it, psmPairs.end()));
 transform(psmPairs.begin(), it, trainPsms.begin(), getPSM);
 transform(it, psmPairs.end(), testPsms.begin(), getPSM);
 if (VERB > 2)
 cerr << "Train set includes now " << trainPsms.size() << " peptides." <<  endl;

 if (!removeDecaying && (VERB > 3))
 {
 cerr << "WARNING: In source fragments will NOT be removed from the test set" << endl;
 cerr << "         Use option -y to remove in source fragments from the test set" << endl;
 }
 else if (VERB > 2)
 cerr << "Test set includes " << testPsms.size() <<  " peptides" << endl;
 }

 // remove non tryptic peptides if the user has wished so
 if (removeNEnzymatic)
 {
 removeNonEnzymatic(trainPsms, "Train");
 if (!testPsms.empty())
 removeNonEnzymatic(testPsms, "Test");
 }

 }*/

void RTPredictor::processTrainData() {
  vector<pair<pair<PSMDescription, bool> , bool> > psmPairs;
  // load input data
  loadTrainFile();
  // remove redundant peptides
  removeRedundantPeptides(trainPsms, "train");
  if (!testFile.empty()) {
    // load test data
    loadTestFile();
    // remove redundant peptides
    if (removeRedundant) {
      removeRedundantPeptides(testPsms, "test");
    } else {
      sort(testPsms.begin(), testPsms.end(), less<PSMDescription> ());
    }
    if (!autoModelSelection) {
      if (VERB > 2) {
        cerr << endl
            << "Removing from the train set the peptides featuring in both train and test sets..."
            << endl;
      }
      vector<PSMDescription>::iterator it;
      int initialSize = trainPsms.size();
      for (it = testPsms.begin(); it != testPsms.end(); ++it) {
        trainPsms.erase(remove(trainPsms.begin(), trainPsms.end(), (*it)),
                        trainPsms.end());
      }
      if (VERB > 2) {
        cerr << (initialSize - trainPsms.size())
            << " peptides were removed." << endl;
        cerr << "Train set includes now " << trainPsms.size()
            << " peptides." << endl;
      }
    }
    // to detect CID fragments
    if (testIncludesRT) {
      addToPairVector(testPsms, false, psmPairs);
    }
  }
  addToPairVector(trainPsms, true, psmPairs);
  removeSourceCIDs(psmPairs);
  // divide in test and train
  if (testPsms.empty() || !testIncludesRT) {
    trainPsms.resize(psmPairs.size());
    transform(psmPairs.begin(), psmPairs.end(), trainPsms.begin(), getPSM);
    if (VERB > 2) {
      cerr << "Train set includes now " << trainPsms.size()
          << " peptides." << endl;
    }
  } else {
    int sizeTest = testPsms.size();
    vector<pair<pair<PSMDescription, bool> , bool> >::iterator it;
    it = partition(psmPairs.begin(), psmPairs.end(), inTrain);
    trainPsms.resize(distance(psmPairs.begin(), it));
    testPsms.resize(distance(it, psmPairs.end()));
    transform(psmPairs.begin(), it, trainPsms.begin(), getPSM);
    transform(it, psmPairs.end(), testPsms.begin(), getPSM);
    if (VERB > 2) {
      cerr << "Train set includes now " << trainPsms.size()
          << " peptides." << endl;
    }
    if (!removeDecaying && (VERB > 3)) {
      cerr
          << "WARNING: In source fragments will NOT be removed from the test set"
          << endl;
      cerr
          << "         Use option -y to remove in source fragments from the test set"
          << endl;
    } else if (VERB > 2) {
      cerr << "Test set includes " << testPsms.size() << " peptides"
          << endl;
    }
  }
  // remove non tryptic peptides if the user has wished so
  if (removeNEnzymatic) {
    removeNonEnzymatic(trainPsms, "train");
    if (!testPsms.empty()) {
      removeNonEnzymatic(testPsms, "test");
    }
  }
}

// memory allocation for the feature table
void RTPredictor::initFeaturesTable(const unsigned int numberRecords,
                                    vector<PSMDescription> & psms,
                                    double* retentionFeatures,
                                    size_t noFeatures) {
  // number of rt related features
  size_t noRTFeatures;
  // copy of the pointer
  double* ptr;
  if (noFeatures == -1) {
    noRTFeatures = model.getNoFeaturesToCalc();
  } else {
    noRTFeatures = noFeatures;
  }
  if (VERB > 3) {
    cerr << endl << "Initializing features table for " << noRTFeatures
        << " features and " << numberRecords << " records..." << endl;
  }
  // allocate memory to all the feature table
  retentionFeatures = new double[noRTFeatures * numberRecords];
  // get the beginning of the memory block
  ptr = retentionFeatures;
  if (!retentionFeatures) {
    if (VERB >= 1) {
      cerr << "Unable to allocate memory for the features table" << endl;
      cerr << "Execution aborted" << endl;
    }
    exit(-1);
  }
  // resize the size of the vector
  psms.resize(numberRecords);
  // divide the allocated memory block to all psms (each one gets a line)
  for (int i = 0; i < numberRecords; i++, ptr += noRTFeatures) {
    psms[i].retentionFeatures = ptr;
  }
  if (VERB > 3) {
    cerr << "Done." << endl;
  }
}

// remove all the duplicate peptides
void RTPredictor::removeRedundantPeptides(vector<PSMDescription> & psms,
                                          const char* dataset) {
  if (VERB > 2) {
    cerr << endl << "Removing redundant peptides from " << dataset
        << "... " << endl;
  }
  int initial, final;
  initial = psms.size();
  // sort first so we keep the peptide with the lowest retention time
  sort(psms.begin(), psms.end(), less<PSMDescription> ());
  psms.resize(distance(psms.begin(), unique(psms.begin(), psms.end())));
  if (VERB > 2) {
    cerr << (initial - psms.size()) << " peptides were removed." << endl;
    cerr << "Dataset includes now " << psms.size() << " peptides." << endl;
  }
}

// remove peptides that are not enzymatic
void RTPredictor::removeNonEnzymatic(vector<PSMDescription> & psms,
                                     string text) {
  vector<PSMDescription>::iterator it1, it2;
  int initialSize;
  initialSize = psms.size();
  if (VERB > 2) {
    cerr << endl << "Remove non enzymatic peptides from " << text << "...";
  }
  it1 = partition(psms.begin(), psms.end(), RTPredictor());
  psms.resize(distance(psms.begin(), it1));
  if (VERB > 2) {
    cerr << (initialSize - psms.size()) << " non-enzymatic peptides "
        << endl;
    cerr << text << " set includes " << psms.size() << " peptides" << endl;
  }
}

// check if a peptides is in source CID (it is included in another peptide, they have similar rt and either
// the father is enzymatic and the child not, either the difference in hydrophobicity is greater than 10
bool RTPredictor::isChildOf(PSMDescription& child, PSMDescription& parent) {
  string peptideP, msPeptideP, peptideC, msPeptideC;
  double rtP, rtC, diff;
  peptideP = parent.getPeptide();
  peptideC = child.getPeptide();
  msPeptideP = getMSPeptide(peptideP);
  msPeptideC = getMSPeptide(peptideC);
  // check sequence inclusion
  if ((msPeptideC.length() >= msPeptideP.length())
      || (msPeptideP.find(msPeptideC) == string::npos)) {
    return false;
  }
  // if the difference in observed rt is greater than 5% of the parent's rt, then the peptide is not aberrant
  // NOT necessary because of the way we search for these peptides
  //rtP = parent.getRetentionTime();
  //rtC = child.getRetentionTime();
  //if ( (abs(rtP - rtC)/rtP) > 0.05 )
  // return false;
  // if the parent is enzymatic and the child is not, then the child is in source CID fragment
  if (theEnzyme->isEnzymatic(peptideP)
      && (!theEnzyme->isEnzymatic(peptideC))) {
    return true;
  }
  // check the difference in hydrophobicities
  diff = RTModel::calcDiffHydrophobicities(msPeptideP, msPeptideC);
  if (diff >= diff_hydrophobicity) {
    return true;
  } else {
    return false;
  }
}

void RTPredictor::addToPairVector(vector<PSMDescription> psms, bool value,
                                  vector<pair<pair<PSMDescription, bool> ,
                                      bool> > & psmPairs) {
  pair<PSMDescription, bool> tmp1;
  pair<pair<PSMDescription, bool> , bool> tmp2;
  for (int i = 0; i < psms.size(); ++i) {
    tmp1.first = psms[i];
    tmp1.second = value;
    tmp2.first = tmp1;
    tmp2.second = false;
    psmPairs.push_back(tmp2);
  }
}

// remove in source CID peptides
void RTPredictor::removeSourceCIDs(vector<pair<
    pair<PSMDescription, bool> , bool> > & psms) {
  double rt, rtp;
  int noPsms = psms.size(), i, j;
  bool isDecay;
  vector<pair<pair<PSMDescription, bool> , bool> >::iterator it;
  if (VERB > 2) {
    if (removeDecaying) {
      cerr << endl
          << "Removing in source fragments from train and test sets..."
          << endl;
    } else {
      cerr << endl << "Removing in source fragments from train set..."
          << endl;
    }
  }
  sort(psms.begin(), psms.end(), mypair);
  for (i = 0; i < noPsms; ++i) {
    rt = psms[i].first.first.getRetentionTime();
    j = i - 1;
    isDecay = false;
    while ((j >= 0) && (((rtp = psms[j].first.first.getRetentionTime())
        * 1.05) >= rt)) {
      if (isChildOf(psms[i].first.first, psms[j].first.first)) {
        psms[i].second = true;
        isDecay = true;
        break;
      }
      j--;
    }
    if (!isDecay) {
      j = i + 1;
      while ((j < noPsms) && (((rtp
          = psms[j].first.first.getRetentionTime()) * 0.95) <= rt)) {
        if (isChildOf(psms[i].first.first, psms[j].first.first)) {
          psms[i].second = true;
          break;
        }
        j++;
      }
    }
  }
  if (!decayingPeptidesFile.empty()) {
    writeDecayToFile(psms);
  }
  int counts = (int)count_if(psms.begin(), psms.end(), decay);
  if (VERB > 2) {
    cerr << counts << " in source fragments were identified." << endl;
  }
  if (removeDecaying) {
    it = remove_if(psms.begin(), psms.end(), decay);
  } else {
    it = remove_if(psms.begin(), psms.end(), decayTrain);
  }
  psms.resize(distance(psms.begin(), it));
}

// write in source CID fragments to a file
void RTPredictor::writeDecayToFile(vector<pair<
    pair<PSMDescription, bool> , bool> > & psms) {
  ofstream outDecay;
  PSMDescription psm;
  outDecay.open(decayingPeptidesFile.c_str());
  outDecay
      << "# File generated by Elude; all decaying peptides are listed below "
      << endl;
  outDecay << "Peptide\tobserved_retention_time\tSet" << endl;
  for (int i = 0; i < psms.size(); ++i)
    if (psms[i].second) {
      psm = psms[i].first.first;
      outDecay << psm.getPeptide() << "\t" << psm.getRetentionTime()
          << "\t";
      if (psms[i].first.second) {
        outDecay << "Train\n";
      } else {
        outDecay << "Test\n";
      }
    }
  outDecay.close();
}

// get the peptide information (for A.MCD.E -> return MCD)
string RTPredictor::getMSPeptide(string& peptide) {
  int pos1, pos2;
  pos1 = peptide.find('.');
  pos2 = peptide.find('.', ++pos1);
  return peptide.substr(pos1, pos2 - pos1);
}

/*
 * Train or load a Support Vector Regressor
 */
void RTPredictor::trainSVRRegressor() {
  double trainingTime = 0.0;
  // memory allocation for the feature table
  initFeaturesTable(trainPsms.size(), trainPsms, trainRetentionFeatures);
  // calculate the hydrophobicity index using the train peptides
  if (model.calculateIndex()) {
    model.computeHydrophobicityIndex(trainPsms);
  }
  // write the generated index to a file
  if (!hydrophobicityFile.empty()) {
    if (model.calculateIndex()) {
      model.printInhouseIndex(hydrophobicityFile);
    } else if (VERB >= 2) {
      cerr << endl
          << "The option too generate a hydrophobicity index is not on, therefore no index is generated."
          << endl;
    }
  }
  // compute the retention features
  model.calcRetentionFeatures(trainPsms);
  // scaling
  // scale both retention times and the values of each retention feature
  PSMDescription::setPSMSet(trainPsms);
  PSMDescription::normalizeRetentionTimes(trainPsms);
  vector<double*> tmp;
  vector<double*> tRetFeat = PSMDescription::getRetFeatures(trainPsms);
  theNormalizer->setSet(tmp,
                        tRetFeat,
                        (size_t)0,
                        model.getNoFeaturesToCalc());
  theNormalizer->normalizeSet(tmp, tRetFeat);
  // train the SVR
  //   clock_t start = std::clock();
  model.trainSVM(trainPsms);
  //   trainingTime = ((std::clock() - start) / (double)CLOCKS_PER_SEC);
  //   if (VERB > 3)
  //        cerr << "Training lasted: "  << trainingTime << " seconds " << endl << endl;
}

// load a SVR Model from a file
void RTPredictor::loadSVRModel() {
  // if a training file was provided and a model already generated, then this step is skipped
  if (trainFile.empty() || model.isModelNull()) {
    if (VERB >= 2) {
      cerr << "Loading model from file " << loadModelFile << endl;
    }
    if (fileExists(loadModelFile)) {
      model.loadSVRModel(loadModelFile, theNormalizer);
      if (VERB >= 2) {
        cerr << "Done" << endl;
      }
    } else {
      if (VERB > 2) {
        cerr << endl << "Error: model file " << loadModelFile
            << " is inaccessible. Please check the file location and permissions and then try again. ";
      }
      exit(1);
    }
  } else {
    if (VERB >= 2) cerr
        << "WARNING: both options -m and -t generate a model.\nSince a model was already generated using data from "
        << trainFile
        << " option -m will be ignored \nPlease use option -m alone to load an existent model\n";
  }
}

/*
 * AUTOMATIC MODEL SELECTION
 */
// read all the files from the libPath
vector<string> RTPredictor::getModelFiles() {
  vector<string> modelFiles;
  DIR* dp;
  struct dirent* dirp;
  string name;
  size_t extension;
  if ((dp = opendir(libPath.c_str())) == NULL) {
    if (VERB >= 1) {
      cerr << "Unable to open " << libPath << endl;
      cerr << "Execution aborted. " << endl;
    }
    exit(-1);
  }
  while ((dirp = readdir(dp)) != NULL) {
    name = string(dirp->d_name);
    extension = name.rfind(".model");
    if (extension == (name.length() - 6)) {
      modelFiles.push_back(libPath + name);
    }
  }
  closedir(dp);
  return modelFiles;
}

// load the model that fits best the training data
void RTPredictor::loadBestModel() {
  // variables to store best correlation
  double bestCorr = 0.0, corr;
  // name of the file including the best model
  string bestModelFile;
  vector<string> modelFiles;
  double* retFeatures;
  bool firstModel = false;
  if (VERB > 2) {
    cerr << "Checking available model files..." << endl;
  }
  // get all the model files from the library
  modelFiles = getModelFiles();
  // if no model files in the library, there is nothing to do at this step
  if (modelFiles.empty()) {
    if (VERB > 1) {
      cerr << "No model file available. " << endl;
    }
    return;
  } else if (VERB > 2) {
    cerr << "The library includes " << modelFiles.size() << " models."
        << endl;
  }
  // if no training data available, we cannot choose a model
  if (trainPsms.empty()) {
    if (VERB > 2) {
      cerr << "Warning: No training data available to select model. "
          << endl;
      cerr << "         The first available model will be selected - "
          << modelFiles[0] << endl;
    }
    bestModelFile = modelFiles[0];
    firstModel = true;
  } else if (modelFiles.size() == 1) {
    bestModelFile = modelFiles[0];
    firstModel = true;
  } else {
    bestModelFile = "";
    for (int i = 0; i < modelFiles.size(); i++) {
      // load the model and information about the normalizer
      if (VERB > 3) {
        cerr << "Processing " << modelFiles[i] << "..." << endl
            << "-------" << endl;
        cerr << "Loading model " << modelFiles[i] << " ..." << endl;
      }
      model.loadSVRModel(modelFiles[i], theNormalizer);
      if (VERB > 3) {
        cerr << "Done." << endl << endl;
      }
      // if the model was not loaded then just give a warning and go to the next one
      if (model.isModelNull()) {
        if (VERB >= 2) {
          cerr << "Warning: Model file " << modelFiles[i]
              << " was not loaded correctly." << endl;
          cerr << "This model will not be considered." << endl;
        }
      } else {
        // normalize the retention times; the scaling parameters should have already been set when the model was loaded
        PSMDescription::normalizeRetentionTimes(trainPsms);
        // calculate the retention features
        model.calcRetentionFeatures(trainPsms);
        // normalize the retention features; the scaling parameters should have already been set when loading the model
        //cerr << "FIXMET: Replace with new normalizer" << endl;
        vector<double*> tmp;
        vector<double*> tRetFeat =
            PSMDescription::getRetFeatures(trainPsms);
        theNormalizer->normalizeSet(tmp, tRetFeat);
        // estimate the retention time
        estimateRetentionTime(trainPsms);
        // compute correlation
        corr = computeSpearmanCorrelation(trainPsms);
        if (corr > bestCorr) {
          bestCorr = corr;
          bestModelFile = modelFiles[i];
        }
        // unnormalize the retention time
        for (int i = 0; i < trainPsms.size(); i++) {
          trainPsms[i].retentionTime
              = PSMDescription::unnormalize(trainPsms[i].retentionTime);
        }
      }
    }
  }
  if (bestModelFile != "") {
    if (VERB > 2) {
      cerr << "----------------------" << endl;
    }
    if (VERB >= 2) {
      cerr << "Model employed " << bestModelFile << endl;
    }
    model.loadSVRModel(bestModelFile, theNormalizer);
    // normalize the retention times; the scaling parameters should have already been set when the model was loaded
    PSMDescription::normalizeRetentionTimes(trainPsms);
    // calculate the retention features
    model.calcRetentionFeatures(trainPsms);
    // normalize the retention features; the scaling parameters should have already been set when loading the model
    //cerr << "FIXMET: Replace with new normalizer" << endl;
    vector<double*> tmp;
    vector<double*> tRetFeat = PSMDescription::getRetFeatures(trainPsms);
    theNormalizer->normalizeSet(tmp, tRetFeat);
    // estimate the retention time
    estimateRetentionTime(trainPsms);
  } else if (VERB >= 2) {
    cerr << "No model loaded " << endl;
  }
}

// find the best line that fits the data (basically the coefficients a, b)
void RTPredictor::findLeastSquaresSolution(
                                           const vector<PSMDescription> & psms,
                                           double& a, double& b) {
  double avgX, avgY, sumX = 0.0, sumY = 0.0, sumXSquared = 0.0, sumXY =
      0.0;
  double ssxx, ssxy;
  int n = psms.size();
  double x, y;
  for (int i = 0; i < n; ++i) {
    y = psms[i].retentionTime;
    x = psms[i].predictedTime;
    //cout << x << ", " << y << endl;
    sumX += x;
    sumY += y;
    sumXSquared += x * x;
    sumXY += x * y;
  }
  avgX = sumX / (double)n;
  avgY = sumY / (double)n;
  ssxx = sumXSquared - (n * avgX * avgX);
  ssxy = sumXY - (n * avgX * avgY);
  a = ssxy / ssxx;
  b = avgY - (a * avgX);
}

// get retention times as vectors
pair<vector<double> , vector<double> > RTPredictor::getRTs(vector<
    PSMDescription> & psms) {
  vector<double> rts, prts;
  for (int i = 0; i < psms.size(); ++i) {
    rts.push_back(psms[i].getRetentionTime());
    prts.push_back(psms[i].getPredictedRetentionTime());
  }
  return make_pair(prts, rts);
}

// load a model from the library
void RTPredictor::loadLibModel() {
  if (!trainFile.empty()) {
    processTrainData();
    // initialize the table storing the retention features
    initFeaturesTable(trainPsms.size(),
                      trainPsms,
                      trainRetentionFeatures,
                      RTModel::totalNumRTFeatures());
  }
  if (VERB > 2) {
    cerr << endl << "--- Load best model for your data ---" << endl;
  }
  // load the model that fits best the data
  loadBestModel();
  // least trimmed regression
  if (linearAdjust && (trainPsms.size() > 0)) {
    unNormalizeRetentionTimes(trainPsms);
    //findLeastSquaresSolution(trainPsms, a, b);
    pair<vector<double> , vector<double> > rts = getRTs(trainPsms);
    lts.setData(rts.first, rts.second);
    if (VERB > 2) {
      cerr << endl << "Compute linear fitting parameters..." << endl;
    }
    lts.runLTS();
  }
}

// add the current model to the library
void RTPredictor::addModelLibrary() {
  if (model.isModelNull()) {
    cerr << endl << "No model available to add to library." << endl;
  } else {
    string libModelFile;
    if (!trainFile.empty()) {
      size_t found;
      found = trainFile.find_last_of("/");
      if (found != string::npos) {
        libModelFile = trainFile.substr(found + 1, trainFile.length());
      } else {
        found = trainFile.find_last_of("\\");
        if (found != string::npos) {
          libModelFile = trainFile.substr(found + 1, trainFile.length());
        } else {
          libModelFile = trainFile;
        }
      }
    } else {
      time_t currentTime;
      struct tm* timeInfo;
      currentTime = time(NULL);
      timeInfo = localtime(&currentTime);
      libModelFile = asctime(timeInfo);
      libModelFile = libModelFile.substr(4, (libModelFile.length() - 10));
      size_t found = libModelFile.find_first_of(" ");
      while (found != string::npos) {
        libModelFile.replace(found, 1, "_");
        found = libModelFile.find_first_of(" ");
      }
    }
    libModelFile = libPath + libModelFile + ".model";
    if (VERB > 2) {
      cerr << endl << "Add current model to library..." << endl;
      cerr << "Model name: " << libModelFile << "." << endl;
    }
    model.saveSVRModel(libModelFile, theNormalizer);
  }
}

/*
 * Testing the model
 */
// load the test data
void RTPredictor::loadTestFile() {
  // retention time
  double retTime;
  string peptide;
  // opens the file for reading
  ifstream in(testFile.c_str(), ios::in);
  string line;
  char* tokens, *tmp;
  if (VERB > 2) {
    cerr << endl << "Loading file " << testFile << "..." << endl;
  }
  // check if the file was open successfully
  if (!in) {
    if (VERB >= 1) {
      cerr << "Unable to open " << testFile << endl;
      cerr << "Please check the file location and try again " << endl;
      cerr << "Execution aborted " << endl;
    }
    exit(-1);
  }
  if (getline(in, line) && !line.empty()) {
    // read the first line and check if the observed rt is included
    tmp = new char[line.size() + 1];
    strcpy(tmp, line.c_str());
    tokens = strtok(tmp, " ");
    peptide.assign(tokens);
    tokens = strtok(NULL, " ");
    if (tokens != NULL) {
      retTime = atof(tokens);
      testPsms.push_back(PSMDescription(peptide, retTime));
    } else {
      testIncludesRT = false;
      testPsms.push_back(PSMDescription(peptide, 0));
    }
    if (testIncludesRT)
      // read the observed retention time and the peptide sequence
      while (in >> peptide >> retTime) {
        testPsms.push_back(PSMDescription(peptide, retTime));
      }
    else
      // read only the peptide
      while (in >> peptide) {
        testPsms.push_back(PSMDescription(peptide, 0));
      }
  }
  if (VERB > 2) {
    cerr << testPsms.size() << " peptides were loaded." << endl;
  }
}

// load and process test data
void RTPredictor::processTestData() {
  // load test data
  loadTestFile();
  // remove redundant peptides
  if (removeRedundant) {
    removeRedundantPeptides(testPsms, "test");
  }
  // remove CID fragments
  if (testIncludesRT
      && (removeDecaying || (!decayingPeptidesFile.empty()))) {
    vector<pair<pair<PSMDescription, bool> , bool> > psmPairs;
    addToPairVector(testPsms, false, psmPairs);
    removeSourceCIDs(psmPairs);
    testPsms.resize(psmPairs.size());
    transform(psmPairs.begin(), psmPairs.end(), testPsms.begin(), getPSM);
  }
  // remove non tryptic peptides
  if (removeNEnzymatic) {
    removeNonEnzymatic(testPsms, "test");
  }
}

// predict retention times using the SVR
void RTPredictor::predictRTs() {
  // memory allocation for the feature table
  initFeaturesTable(testPsms.size(), testPsms, testRetentionFeatures);
  // compute the retention features
  model.calcRetentionFeatures(testPsms);
  // scaling
  // retention times and all the features are scaled
  PSMDescription::normalizeRetentionTimes(testPsms);
  vector<double*> tmp;
  vector<double*> tRetFeat = PSMDescription::getRetFeatures(testPsms);
  theNormalizer->normalizeSet(tmp, tRetFeat);
  // predict rts
  estimateRetentionTime(testPsms);
  unNormalizeRetentionTimes(testPsms);
}

// estimates the retention time and fills in the predicted retention time field in each PSMDescription
void RTPredictor::estimateRetentionTime(vector<PSMDescription> & psms) {
  assert(model.getModel());
  vector<PSMDescription>::iterator it;
  double predictedRT;
  if (VERB > 2) {
    cerr << endl << "Predicting retention times..." << endl;
  }
  for (it = psms.begin(); it != psms.end(); ++it) {
    predictedRT = model.estimateRT(it->retentionFeatures);
    //cout << predictedRT << " ";
    it->predictedTime = predictedRT;
  }
  if (VERB > 2) {
    cerr << "Done." << endl;
  }
}

// write predictions to the output file
void RTPredictor::writeOutputFile(vector<PSMDescription> & psms) {
  ofstream out(outputFile.c_str());
  vector<PSMDescription>::iterator it;
  double predictedRT, observedRT;
  out << "# File generated by Elude " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  if (testIncludesRT) {
    out << "Peptide\tpredicted_retention_time\t[observed_retention_time]"
        << endl;
  } else {
    out << "Peptide\tpredicted_retention_time" << endl;
  }
  for (it = psms.begin(); it != psms.end(); ++it) {
    if (testIncludesRT) {
      out << (it->peptide) << "\t" << it->getPredictedRetentionTime()
          << "\t" << it->getRetentionTime() << "\n";
    } else {
      out << (it->peptide) << "\t" << it->getPredictedRetentionTime()
          << "\n";
    }
  }
  out.close();
}

/*
 void RTPredictor::run()
 {
 LTSRegression lts;
 vector<double> x;
 vector<double> y;

 for(int i = 0; i < 12000; ++i)
 {
 x.push_back((double)i*100);
 y.push_back((double)i);
 }

 for(int i = 0; i < 8000; ++i)
 {
 x.push_back((double)i);
 y.push_back((double)i*20);
 }
 for(int i = 0; i < 200; ++i)
 {
 x.push_back((double)i);
 y.push_back((double)i*50);
 }
 lts.setData(x, y);
 lts.runLTS();
 pair<double, double> p = lts.getRegCoefficients();
 cout << "y = " << p.first << " * x + " << p.second << endl;
 }*/

/*
 * RUN the application
 */
// run Elude
void RTPredictor::run() {
  // print information about the current run
  if (VERB > 2) {
    printRunInformation();
  }
  // get a normalizer
  theNormalizer = Normalizer::getNormalizer();
  theNormalizer->resizeVecs(RTModel::totalNumRTFeatures());
  if (autoModelSelection)
  // select the model that fits best the training data
  {
    loadLibModel();
  } else {
    if (!trainFile.empty()) {
      // train a model
      processTrainData();
      //writeRetentionTimeFile("../results/article_results/data/clean_train/Luna_60_clean_train.txt" , trainPsms);
      trainSVRRegressor();
    }
    if (!loadModelFile.empty())
    // load a model from a file
    {
      loadSVRModel();
    }
  }
  if (!testFile.empty()) {
    double pcorrelation, scorrelation, ms, timeWindow;
    // check that a retention model was generated
    if (model.isModelNull()) {
      if (VERB >= 1) {
        cerr << "Error: No SVM model available " << endl;
        cerr
            << "Please train a model using the -t option or load an existing model from a file using the"
              " -m option and try again" << endl;
        cerr << "Execution aborted" << endl;
      }
      exit(-1);
    }
    if (trainFile.empty()) {
      // process test data if not already generated
      processTestData();
    }
    // predict retention time(results are unnormalized)
    predictRTs();
    // calibration (only if the model was select automatically from the library)
    if (autoModelSelection && linearAdjust && (trainPsms.size() > 0)) {
      if (VERB > 2) {
        cerr << "Adjusting linearly..." << endl;
      }
      for (int i = 0; i < testPsms.size(); i++) {
        testPsms[i].predictedTime = lts.predict(testPsms[i].predictedTime);
      }
      if (VERB > 2) {
        cerr << "Done." << endl;
      }
    }
    if (testIncludesRT) {
      // compute correlations and ms
      pcorrelation = computePearsonCorrelation(testPsms);
      scorrelation = computeSpearmanCorrelation(testPsms);
      timeWindow = computeWindow(testPsms);
      if (VERB >= 2) {
        cerr << endl << "--------------------------------------" << endl;
        cerr << "Pearson Correlation = " << pcorrelation << endl;
        cerr << "Spearman Correlation = " << scorrelation << endl;
        cerr << fractionPeptides * 100 << "% window = " << timeWindow
            << endl;
        cerr << "--------------------------------------" << endl;
      }
    }
    // write the output file
    if (!outputFile.empty()) {
      writeOutputFile(testPsms);
    } else {
      cerr << endl << "Predictions: " << endl;
      cerr << "Peptide Predicted_retention_time" << endl;
      for (int i = 0; i < testPsms.size(); ++i) {
        cerr << testPsms[i].peptide << " " << testPsms[i].predictedTime
            << endl;
      }
    }
  }
  // save the model
  if (!saveModelFile.empty()) {
    if (model.isModelNull()) {
      cerr
          << "WARNING: No model was generated; please use options -m and -t to generate a model and then save it. \n"
          << "Option -s has no effect;" << endl;
    } else {
      // SAVE A MODEL in the given file
      cerr << endl << "Saving SVR model to file " << saveModelFile << endl;
      model.saveSVRModel(saveModelFile, theNormalizer);
      cerr << "Done." << endl;
    }
  }
  // adsd the current model to the library
  if (addLibModel) {
    addModelLibrary();
  }
}

// unnormalize both observed and predicted retention times for a set of peptides
void RTPredictor::unNormalizeRetentionTimes(vector<PSMDescription> & psms) {
  for (int i = 0; i < psms.size(); ++i) {
    psms[i].retentionTime
        = PSMDescription::unnormalize(psms[i].retentionTime);
    psms[i].predictedTime
        = PSMDescription::unnormalize(psms[i].predictedTime);
  }
}

/*
 * Compute MS and correlations
 */

// calculate ms between observed and predicted retention times
double RTPredictor::computeMS(vector<PSMDescription> & psms) {
  double ms = 0.0, diff = 0.0;
  int noPsms;
  if (VERB > 2) {
    cerr << endl << "Computing ms..." << endl;
  }
  noPsms = psms.size();
  for (int i = 0; i < noPsms; i++) {
    diff = psms[i].predictedTime - psms[i].retentionTime;
    ms += diff * diff;
  }
  ms = ms / noPsms;
  if (VERB > 2) {
    cerr << "Done." << endl;
  }
  return ms;
}

// calculate the correlation rho
double RTPredictor::computePearsonCorrelation(
                                              vector<PSMDescription> & psms) {
  double sumObserved = 0.0, sumPredicted = 0.0;
  double meanObserved = 0.0, meanPredicted = 0.0;
  double stdevObserved = 1, stdevPredicted = 1;
  double devObserved, devPredicted;
  double numerator = 0.0;
  double corr;
  int noPsms = psms.size();
  if (VERB > 2) {
    cerr << endl << "Computing Pearson correlation..." << endl;
  }
  // calculate means
  for (int i = 0; i < noPsms; i++) {
    sumObserved += psms[i].retentionTime;
    sumPredicted += psms[i].predictedTime;
  }
  meanObserved = sumObserved / noPsms;
  meanPredicted = sumPredicted / noPsms;
  sumObserved = 0.0;
  sumPredicted = 0.0;
  // calculate stdevs
  for (int i = 0; i < noPsms; i++) {
    devObserved = psms[i].retentionTime - meanObserved;
    devPredicted = psms[i].predictedTime - meanPredicted;
    numerator += devObserved * devPredicted;
    sumObserved += pow(devObserved, 2);
    sumPredicted += pow(devPredicted, 2);
  }
  stdevObserved = sqrt(sumObserved / (noPsms - 1));
  stdevPredicted = sqrt(sumPredicted / (noPsms - 1));
  corr = numerator / ((noPsms - 1) * stdevObserved * stdevPredicted);
  if (VERB > 2) {
    cerr << "r_pearson = " << corr << endl;
    cerr << "Done." << endl;
  }
  return corr;
}

// calculate the spearman correlation
double RTPredictor::computeSpearmanCorrelation(
                                               vector<PSMDescription> & psms) {
  double corr = 0.0, d = 0.0, avgRank, rankP, rankO;
  int i, j;
  int n = psms.size();
  vector<pair<PSMDescription, double> > rankedPsms;
  if (VERB > 2) {
    cerr << endl << "Computing Spearman correlation..." << endl;
  }
  // sort peptides according to observed retention time
  sort(psms.begin(), psms.end(), comparePsmsRT);
  // record ranks
  i = 0;
  while (i < n) {
    avgRank = j = i + 1;
    while ((j < n) && (psms[i].getRetentionTime()
        == psms[j].getRetentionTime())) {
      avgRank += ++j;
    }
    avgRank = avgRank / (double)(j - i);
    for (int k = i; k < j; ++k) {
      rankedPsms.push_back(make_pair(psms[k], avgRank));
    }
    i = j;
  }
  // sort peptides according to predicted rt
  sort(rankedPsms.begin(), rankedPsms.end(), comparePsmsPRT);
  // calculate sum of squared differences btw ranks
  i = 0;
  while (i < n) {
    // calculate rank of predicted rt
    rankP = j = i + 1;
    while ((j < n) && (rankedPsms[i].first.getPredictedRetentionTime()
        == rankedPsms[j].first.getPredictedRetentionTime())) {
      rankP += ++j;
    }
    rankP = rankP / (double)(j - i);
    // calculate and add squared difference
    for (int k = i; k < j; ++k) {
      d += pow(rankedPsms[k].second - rankP, 2);
    }
    // increase i
    i = j;
  }
  corr = 1.0 - ((6.0 * d) / (double)(n * (pow(n, 2.) - 1)));
  if (VERB > 2) {
    cerr << "r_spearman = " << corr << endl;
    cerr << "Done." << endl;
  }
  return corr;
}

/*
 * Calculate a window where fraction of the peptides are included
 */

/*
 double RTPredictor::computeWindow(vector<PSMDescription> & psms)
 {
 double win, diff;
 int nr = (int)round(fractionPeptides * (double)psms.size());
 int k = psms.size() - nr, i = 1;

 sort(psms.begin(), psms.end(), comparePsmsDeltaRT);

 win = (PSMDescription::unnormalize(psms[nr - 1].getPredictedRetentionTime()) - PSMDescription::unnormalize(psms[nr - 1].getRetentionTime())) -
 (PSMDescription::unnormalize(psms[0].getPredictedRetentionTime()) - PSMDescription::unnormalize(psms[0].getRetentionTime()));

 while (i <= k)
 {
 diff = (PSMDescription::unnormalize(psms[i + nr - 1].getPredictedRetentionTime()) - PSMDescription::unnormalize(psms[i + nr - 1].getRetentionTime())) -
 (PSMDescription::unnormalize(psms[i].getPredictedRetentionTime()) - PSMDescription::unnormalize(psms[i].getRetentionTime()));
 if (diff < win)
 win = diff;
 i++;
 }

 return win;
 }*/

double RTPredictor::computeWindow(vector<PSMDescription> & psms) {
  double win, diff;
  int nr = (int)round(fractionPeptides * (double)psms.size());
  int k = psms.size() - nr, i = 1;
  if (VERB > 2) {
    cerr << endl << "Computing " << fractionPeptides * 100
        << "% window..." << endl;
  }
  sort(psms.begin(), psms.end(), comparePsmsDeltaRT);
  win = (psms[nr - 1].getPredictedRetentionTime()
      - psms[nr - 1].getRetentionTime())
      - (psms[0].getPredictedRetentionTime() - psms[0].getRetentionTime());
  while (i <= k) {
    diff = (psms[i + nr - 1].getPredictedRetentionTime()
        - psms[i + nr - 1].getRetentionTime())
        - (psms[i].getPredictedRetentionTime()
            - psms[i].getRetentionTime());
    if (diff < win) {
      win = diff;
    }
    i++;
  }
  if (VERB > 2) {
    cerr << "Done." << endl;
  }
  return win;
}

bool RTPredictor::fileExists(string& strFilename) {
  struct stat stFileInfo;
  int intStat;
  if (stat(strFilename.c_str(), &stFileInfo) != 0) {
    return false;
  } else {
    return true;
  }
}
/*
 * Printing functions
 */
// write peptides and retention times to a file
void RTPredictor::writeRetentionTimeFile(const char* fileName, vector<
    PSMDescription> & psms) {
  ofstream out(fileName);
  vector<PSMDescription>::iterator it;
  for (it = psms.begin(); it != psms.end(); ++it) {
    out << (it->peptide) << " " << (it->retentionTime) << "\n";
  }
  out.close();
}

// write the value of the features to a file
void RTPredictor::writeFeaturesFile(const char* file, vector<
    PSMDescription> & psms, bool unnormalized) {
  ofstream out(file);
  vector<PSMDescription>::iterator it;
  double predictedRT, observedRT;
  out << "# File generated by Elude " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  out << "Peptide\tpredicted_retention_time\tobserved_retention_time"
      << endl;
  if (unnormalized) {
    vector<double*> tmp = PSMDescription::getRetFeatures(psms);
    theNormalizer->unNormalizeSet(tmp);
  }
  for (it = psms.begin(); it != psms.end(); ++it) {
    predictedRT = it->predictedTime;
    observedRT = it->retentionTime;
    out << (it->peptide) << "\t" << predictedRT << "\t" << observedRT;
    for (size_t i = 0; i < *(theNormalizer->getNumRetFeatures()); ++i) {
      out << "\t" << it->retentionFeatures[i];
    }
    out << "\n";
  }
  out.close();
}

// print a vector of PSMs
void RTPredictor::printPsms(vector<PSMDescription> & psms) {
  cerr << endl << "Printing..." << endl;
  vector<PSMDescription>::iterator it;
  for (it = psms.begin(); it != psms.end(); it++) {
    cerr << it->peptide << it->retentionTime << endl;
  }
  cerr << "Done." << endl << endl;
}
