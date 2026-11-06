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
 * direction is not worse than the initial one: by the targets with q < test
 * FDR summed over the test bins, each with its own q values, or by those of
 * the merged out-of-fold list, see DirectionCheck and Caller's
 * --direction-check.
 */
class SanityCheck {
 public:
  // How validateDirection counts the targets with q < test FDR. With one
  // bin, as with RESET, all three are the same.
  enum class DirectionCheck {
    BINS,           // within each test bin, with the bin's own q values
    SECOND_CHANCE,  // as BINS; where that count would reset, the trained
                    // SVMs are kept if the merged out-of-fold list,
                    // rescaled by the training decoys, holds targets with
                    // q < test FDR and no fewer than for the initial
                    // direction counted in the same way
    POOLED          // in the merged out-of-fold list only, merged with the
                    // scheme of setMergeScheme
  };

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
  void setDirectionCheck(DirectionCheck check) { directionCheck_ = check; }
  DirectionCheck directionCheck() const { return directionCheck_; }
  // whether the last validateDirection kept the trained SVMs by the second
  // chance, after the per-bin count would have reset them
  bool secondChanceKept() const { return secondChanceKept_; }
  // the scheme of the final merge when the second chance keeps the trained
  // SVMs; train-decoys, as counted, unless the caller sets another
  void setSecondChanceMergeScheme(MergeNormalization scheme) {
    secondChanceMergeScheme_ = scheme;
  }
  MergeNormalization secondChanceMergeScheme() const {
    return secondChanceMergeScheme_;
  }
  // how the test bins are merged, which the pooled check repeats, see
  // Scores::merge and Scores::countMergedTargets; the second chance takes
  // fdr and skipNormalizeScores from here and always counts with
  // train-decoys
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
    return directionCheck_ == DirectionCheck::POOLED && numBins > 1u;
  }
  // whether a failed per-bin count of these bins gets the second chance
  bool useSecondChance(std::size_t numBins) const {
    return directionCheck_ == DirectionCheck::SECOND_CHANCE && numBins > 1u;
  }
  // whether the weight of the --default-direction feature has the wrong
  // sign in one of the folds' vectors w
  bool wrongSignOfDefaultDirection(
      const std::vector<std::vector<double> >& w) const;
  // The second chance: the targets with q < test FDR in the merged
  // out-of-fold list with train-decoys for the trained weights w and for
  // the initial direction of getInitDirection, and whether they keep w,
  // which also needs the right sign of the --default-direction feature.
  // Prints its decision.
  bool secondChance(const std::vector<Scores>& testset,
                    const std::vector<Scores>& trainset,
                    const std::vector<std::vector<double> >& w,
                    int& trainedCount) const;

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
  std::vector<double> fileWeights_;
  // input from concatenated search, i.e. not separate target and decoy searches
  bool concatenatedSearch_ = true;
  // the direction check and the merge that the pooled count repeats
  DirectionCheck directionCheck_ = DirectionCheck::BINS;
  MergeNormalization mergeScheme_ = MergeNormalization::TEST;
  double mergeFdr_ = 0.01;
  bool skipNormalizeScores_ = false;
  // the second chance: the initial direction of each fold, its decision in
  // the last validateDirection and the scheme of the merge it implies
  std::vector<std::vector<double> > initWeights_;
  bool secondChanceKept_ = false;
  MergeNormalization secondChanceMergeScheme_ =
      MergeNormalization::TRAIN_DECOYS;
};

#endif /*SANITYCHECK_H_*/
