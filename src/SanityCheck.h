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
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_

#include <iostream>
#include <string>
#include <vector>

#include "Scores.h"

class Normalizer;

/*
 * SanityCheck sets the initial direction of the SVM, either from an input
 * weight file, a default direction given by the user or in the input file,
 * or the single most discriminative feature, and checks that the trained
 * direction is not worse than the initial one: either by the targets with
 * q < test FDR summed over the test bins, each with its own q values (bins,
 * the default), or by those of the merged out-of-fold list (pooled, see
 * Caller's --direction-check).
 */
class SanityCheck {
 public:
  void readWeights(std::istream& weightStream, std::vector<double>& w);
  void checkAndSetDefaultDir();
  int getInitDirection(std::vector<Scores>& testset,
                       std::vector<Scores>& trainset, const Normalizer* pNorm,
                       std::vector<std::vector<double> >& w, double test_fdr,
                       double initial_train_fdr);
  bool validateDirection(std::vector<Scores>& testset,
                         std::vector<Scores>& trainset,
                         std::vector<std::vector<double> >& w);
  void resetDirection(std::vector<Scores>& trainset,
                      std::vector<std::vector<double> >& w);

  void setInitWeightFN(const std::string& fn) { initWeightFN_ = fn; }
  void setInitDefaultDirName(const std::string& dirName) {
    initDefaultDirName_ = dirName;
  }
  void setOverrule(bool orl) { overRule_ = orl; }
  // with --override a failed check keeps the trained direction
  bool getOverrule() const { return overRule_; }
  // count the targets with q < test FDR in the merged out-of-fold list
  // instead of in each test bin; with one bin both are the same
  void setPooledCheck(bool pooled) { pooledCheck_ = pooled; }
  bool pooledCheck() const { return pooledCheck_; }
  // how the test bins are merged, which the pooled check repeats, see
  // Scores::merge and Scores::countMergedTargets
  void setMergeScheme(MergeNormalization scheme, double fdr,
                      bool skipNormalizeScores) {
    mergeScheme_ = scheme;
    mergeFdr_ = fdr;
    skipNormalizeScores_ = skipNormalizeScores;
  }
  void addDefaultWeights(const std::vector<double>& defaultWeights) {
    defaultWeights_ = defaultWeights;
  }
  const std::vector<double>& getDefaultWeights() const {
    return defaultWeights_;
  }
  bool concatenatedSearch() const { return concatenatedSearch_; }
  void setConcatenatedSearch(const bool s) { concatenatedSearch_ = s; }

 protected:
  void getDefaultDirection(std::vector<Scores>& trainset,
                           std::vector<std::vector<double> >& w);
  // whether the targets are counted in the merged list of these bins
  bool usePooledCount(std::size_t numBins) const {
    return pooledCheck_ && numBins > 1u;
  }

  int initPositives_ = 0;
  double test_fdr_ = 0.01, initial_train_fdr_ = 0.01;
  bool overRule_ = false;
  std::string initWeightFN_;
  std::string initDefaultDirName_;
  // Default Direction, 0=do not use, positive integer = feature number,
  // negative integer = lower score better
  int initDefaultDir_ = 0;
  // initial direction given in the input file
  std::vector<double> defaultWeights_;
  // normalized initial weights read from the weights file (-W), if any
  std::vector<double> initWeights_;
  // input from concatenated search, i.e. not separate target and decoy searches
  bool concatenatedSearch_ = true;
  // the pooled direction check and the merge it repeats
  bool pooledCheck_ = false;
  MergeNormalization mergeScheme_ = MergeNormalization::TEST;
  double mergeFdr_ = 0.01;
  bool skipNormalizeScores_ = false;
};

#endif /*SANITYCHECK_H_*/
