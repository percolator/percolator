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
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"

void SanityCheck::checkAndSetDefaultDir() {
  if (!initDefaultDir_ && initDefaultDirName_.size() > 0) {
    int sign = 1;
    if (initDefaultDirName_[0] == '-') {
      initDefaultDirName_.erase(0,1);
      sign = -1;
    }
    initDefaultDir_ = sign * DataSet::getFeatureNames().getFeatureNumber(initDefaultDirName_);
    if (initDefaultDir_ == 0) {
      ostringstream temp;
      temp << "ERROR: Initial direction feature name \"" << initDefaultDirName_ << "\" not found" << std::endl;
      throw MyException(temp.str());
    }
  }
}

int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  const Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr,
                                  double initial_train_fdr) {
  test_fdr_= test_fdr;
  initial_train_fdr_ = initial_train_fdr;
  if (initWeightFN_.size() > 0) {
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN_.data(), ios::in);
    if (weightStream.is_open()) {
      readWeights(weightStream, ww);
      weightStream.close();
      assert(pNorm);
      pNorm->normalizeweight(ww, w[0]);
      for (size_t set = 1; set < w.size(); ++set) {
        copy(w[0].begin(), w[0].end(), w[set].begin());
      }
      fileWeights_ = w[0];
    } else {
      std::cerr << "WARNING: Could not find weights input file " << initWeightFN_
                << ". Using default weights instead." << std::endl;
      getDefaultDirection(trainset, w);
    }
  } else {
    getDefaultDirection(trainset, w);
  }
  initPositives_ = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    initPositives_ += testset[set].calcScoresAndQvals(w[set], test_fdr);
  }
  // The pooled check compares the trained direction with the initial one
  // counted in the same way: per fold, merged and pooled. The bins are
  // still scored above, so that both checks leave them in the same state.
  if (usePooledCount(w.size())) {
    initPositives_ = std::max(
        0, Scores::countMergedTargets(testset, trainset, w, mergeScheme_,
                                      mergeFdr_, skipNormalizeScores_,
                                      test_fdr));
  }
  // The second chance counts the initial direction in the merged list only
  // when it runs; training changes w, so the vectors are kept for it.
  initWeights_.clear();
  if (useSecondChance(w.size())) {
    initWeights_ = w;
  }
  return initPositives_;
}

void SanityCheck::getDefaultDirection(vector<Scores>& trainset,
                                      vector<vector<double> >& w) {
    
  //If I have not been given a initial direction
  if (!initDefaultDir_) {
    //if the default weights from the input file are not present
    if (defaultWeights_.size() == 0) {
      // Set init direction to be the most discriminative direction
      for (size_t set = 0; set < w.size(); ++set) {
        if (VERB > 1) std::cerr << "Split " << set + 1 << ":\t";
        trainset[set].getInitDirection(initial_train_fdr_, w[set]);
      }
    } else {
      // I want to assign the default vector that is present in the input file
      for (size_t set = 0; set < w.size(); ++set) {
        for (size_t ix = 0; ix < w[set].size(); ix++) {
          w[set][ix] = 0;
          if (ix < defaultWeights_.size()){
              w[set][ix] = defaultWeights_[ix];
          }          
        }
      }
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      for (size_t ix = 0; ix < w[set].size(); ix++) {
        w[set][ix] = 0;
      }
      w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] = (initDefaultDir_ < 0 ? -1 : 1);
    }
  }
}

bool SanityCheck::wrongSignOfDefaultDirection(
    const vector<vector<double> >& w) const {
  if (!initDefaultDir_) {
    return false;
  }
  for (size_t set = 0; set < w.size(); ++set) {
    if (w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] *
            initDefaultDir_ <= 0) {
      return true;
    }
  }
  return false;
}

bool SanityCheck::secondChance(const vector<Scores>& testset,
                               const vector<Scores>& trainset,
                               const vector<vector<double> >& w,
                               int& trainedCount) const {
  // the anchors of the training sets exist for test bins of any size
  const MergeNormalization counted = MergeNormalization::TRAIN_DECOYS;
  trainedCount =
      Scores::countMergedTargets(testset, trainset, w, counted, mergeFdr_,
                                 skipNormalizeScores_, test_fdr_);
  int initCount = 0;
  if (initWeights_.size() == w.size()) {
    initCount = std::max(
        0, Scores::countMergedTargets(testset, trainset, initWeights_,
                                      counted, mergeFdr_,
                                      skipNormalizeScores_, test_fdr_));
  }
  const bool countsKeep = trainedCount > 0 && trainedCount >= initCount;
  // the sign check of --default-direction applies as after a passed count
  const bool wrongSign = countsKeep && wrongSignOfDefaultDirection(w);
  const bool keep = countsKeep && !wrongSign;
  cerr << "Second chance on the merged out-of-fold list"
       << (skipNormalizeScores_ ? "" : ", rescaled by the training decoys")
       << ": ";
  if (trainedCount < 0) {
    cerr << "the scores of a test bin cannot be rescaled with the trained "
         << "weights";
  } else {
    cerr << trainedCount << " targets with q<" << test_fdr_
         << " for the trained SVMs, " << initCount
         << " for the initial direction";
  }
  if (overRule_) {
    cerr << "; --override keeps the trained SVMs in any case." << endl;
  } else if (keep) {
    cerr << "; the trained SVMs are kept";
    if (!skipNormalizeScores_) {
      cerr << ", merged with "
           << mergeNormalizationName(secondChanceMergeScheme_);
    }
    cerr << "." << endl;
  } else if (wrongSign) {
    cerr << "; the trained SVMs are not kept, as the weight of the main "
         << "scoring direction (--default-direction) has the wrong sign."
         << endl;
  } else {
    cerr << "; the trained SVMs are not kept." << endl;
  }
  return keep;
}

bool SanityCheck::validateDirection(vector<Scores>& testset,
                                    vector<Scores>& trainset,
                                    vector<vector<double> >& w) {
  secondChanceKept_ = false;
  const bool pooled = usePooledCount(w.size());
  int overFDR = 0;
  if (pooled) {
    // on copies: the test bins are only scored once the direction is decided
    overFDR = Scores::countMergedTargets(testset, trainset, w, mergeScheme_,
                                         mergeFdr_, skipNormalizeScores_,
                                         test_fdr_);
    if (overFDR < 0) {
      cerr << "The scores of a test bin cannot be rescaled to merge the "
           << "bins with the trained weights" << endl;
      resetDirection(trainset, w);
      return false;
    }
    if (VERB > 1) {
      cerr << "Found " << overFDR << " targets with q<" << test_fdr_
           << " in the merged out-of-fold list, " << initPositives_
           << " in the initial direction" << endl;
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      overFDR += testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  bool keptBySecondChance = false;
  if (overFDR <= 0 || initPositives_ > overFDR) {
    if (overFDR <= 0) {
      cerr << "No targets found with q<" << test_fdr_ << endl;
    } else {
      cerr << "Less identifications (" << overFDR << " vs " << initPositives_
          << ") after percolator processing than before processing" << endl;
    }
    // Where the per-bin count would reset, the second chance judges the
    // trained SVMs on the merged out-of-fold list, which small test bins
    // do not limit. It decides nothing under --override.
    const bool secondChanceTried = useSecondChance(w.size());
    int pooledCount = 0;
    keptBySecondChance = secondChanceTried &&
                         secondChance(testset, trainset, w, pooledCount) &&
                         !overRule_;
    if (!keptBySecondChance) {
      resetDirection(trainset, w);
      if (overFDR <= 0 && !overRule_) {
        if (pooled || (secondChanceTried && pooledCount <= 0)) {
          cerr << "Hint: q<" << test_fdr_ << " needs about "
               << 1.0 / test_fdr_ << " targets above the best-scoring decoy "
               << "of the merged list; a larger --testFDR is an alternative."
               << endl;
        } else if (w.size() > 1u && !secondChanceTried) {
          cerr << "Hint: each of the " << w.size() << " test folds may be too "
               << "small to reach q<" << test_fdr_ << ". --num-folds 1 trains "
               << "and evaluates on the same PSMs," << endl
               << "so its q values and PEPs are optimistic, and wrong for "
               << "inputs without signal; a larger --testFDR is an "
               << "alternative." << endl;
        }
      }
      return false;
    }
  }
  if (wrongSignOfDefaultDirection(w)) {
    if (VERB > 1) {
      cerr << "Warning, wrong sign of the weight for main scoring direction"
           << endl;
    }
    resetDirection(trainset, w);
    return false;
  }
  if (pooled) {
    // the direction is kept: score the test bins with their folds' weights
    // for the merge, as the count of the bins check does
    for (size_t set = 0; set < w.size(); ++set) {
      testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  // the per-bin count has already scored the test bins with w
  secondChanceKept_ = keptBySecondChance;
  return true;
}

/**
 * Reads the raw weights of the first split from a file written with -w:
 * comment lines, then per split a line of feature names followed by m0, a
 * line of normalized weights and a line of raw weights. The weights are
 * used by position, so the names, if present, must be the feature names of
 * the input, and the raw weights line must hold one finite number per
 * feature plus m0.
 */
void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  const size_t numWeights = FeatureNames::getNumFeatures() + 1;
  std::string line, namesLine;
  unsigned int lineNr = 0;
  bool foundWeights = false;
  // Find first line with weights, but skip this and get second line containing raw features
  while (std::getline(weightStream, line)) {
    ++lineNr;
    if (line[0] == '-' || (line[0] >= '0' && line[0] <= '9')) {
      foundWeights = true;
      break;
    }
    if (!line.empty() && line[0] != '#') namesLine = line;
  }
  std::string rawLine;
  if (!foundWeights || !std::getline(weightStream, rawLine)) {
    throw MyException("ERROR: The weights file " + initWeightFN_ +
                      " has no line of normalized weights followed by a "
                      "line of raw weights.\n");
  }
  if (!namesLine.empty()) {
    std::vector<std::string> names;
    std::istringstream nameStream(namesLine);
    for (std::string name; std::getline(nameStream, name, '\t');) {
      names.push_back(name);
    }
    while (!names.empty() && names.back().find_first_not_of(" \r") ==
                                 std::string::npos) {
      names.pop_back();  // a trailing tab or carriage return
    }
    if (!names.empty()) {
      size_t ix = names.back().find_last_not_of(" \r");
      names.back().erase(ix + 1);
    }
    // the names as -w writes them for this input
    std::vector<std::string> expected;
    std::istringstream expectedStream(
        DataSet::getFeatureNames().getFeatureNames() + "\tm0");
    for (std::string name; std::getline(expectedStream, name, '\t');) {
      expected.push_back(name);
    }
    for (size_t ix = 0; ix < std::max(names.size(), expected.size()); ++ix) {
      if (ix < names.size() && ix < expected.size() &&
          names[ix] == expected[ix]) {
        continue;
      }
      ostringstream oss;
      oss << "ERROR: The weights file " << initWeightFN_ << " has "
          << names.size() << " columns, but the input has "
          << expected.size() << " (its features and m0). Column " << ix + 1
          << " is " << (ix < names.size() ? "\"" + names[ix] + "\"" :
                                            std::string("missing"))
          << " in the weights file and "
          << (ix < expected.size() ? "\"" + expected[ix] + "\"" :
                                     std::string("missing"))
          << " in the input. Weights are used by the position of the "
          << "feature, so the file must come from an input with the same "
          << "features in the same order." << std::endl;
      throw MyException(oss.str());
    }
  }
  std::istringstream rawStream(rawLine);
  std::vector<double> raw;
  for (double value; rawStream >> value;) raw.push_back(value);
  bool allFinite = std::all_of(raw.begin(), raw.end(),
                               [](double x) { return isfinite(x); });
  if (!rawStream.eof() || raw.size() != numWeights || !allFinite) {
    ostringstream oss;
    oss << "ERROR: Line " << lineNr + 1 << " of the weights file "
        << initWeightFN_ << " should hold " << numWeights
        << " finite raw weights, one per feature and m0, but ";
    if (!rawStream.eof() || !allFinite) {
      oss << "holds a value that is not a finite number";
    } else {
      oss << "holds " << raw.size();
    }
    oss << "." << std::endl;
    throw MyException(oss.str());
  }
  std::copy(raw.begin(), raw.end(), w.begin());
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<Scores>& trainset,
                                 vector<vector<double> >& w) {
  if (!overRule_ && !fileWeights_.empty()) {
    // the vector given with -W, not the default that it replaced
    cerr << "Resetting score vector, using the initial weights from the "
         << "weights file. Use --override flag to prevent this." << endl;
    for (size_t set = 0; set < w.size(); ++set) {
      w[set] = fileWeights_;
    }
  } else if (!overRule_) {
    cerr << "Resetting score vector, using default vector. Use --override flag to prevent this." << endl;
    getDefaultDirection(trainset, w);
  }
}

