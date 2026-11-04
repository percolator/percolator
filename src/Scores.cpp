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

#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "DataSet.h"
#include "Globals.h"
#include "MassHandler.h"
#include "Normalizer.h"
#include "PosteriorEstimator.h"
#include "Scores.h"
#include "SetHandler.h"
#include "SequentialSort.h"
#include "ssl.h"
#include "IsotonicPEP.h"

namespace {
inline bool inParallelRegion() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}
}  // namespace

void Scores::merge(std::vector<Scores>& sv,
                   double fdr,
                   bool skipNormalizeScores,
                   std::vector<std::vector<double> >& all_w) {
  scores_.clear();
  size_t numScores = 0;
  for (const Scores& binScores : sv) {
    numScores += binScores.scores_.size();
  }
  scores_.reserve(numScores);
  std::vector<Scores>::iterator cvBinScores = sv.begin();
  std::vector<std::vector<double> >::iterator weights = all_w.begin();
  for (; cvBinScores != sv.end(); cvBinScores++, weights++) {
    // the test sets are usually sorted already, by calcScores
    if (!cvBinScores->isSortedByScore()) {
      sort(cvBinScores->begin(), cvBinScores->end(), greater<ScoreHolder>());
    }
    cvBinScores->checkSeparationAndSetPi0();
    cvBinScores->calcQvals(fdr);
    if (!skipNormalizeScores) {
      cvBinScores->normalizeScores(fdr, *weights, sv.size() == 1u);
    }
    copy(cvBinScores->begin(), cvBinScores->end(), back_inserter(scores_));
    // release the bin before the sort of all PSMs
    cvBinScores->reset();
  }
  postMergeStep();
}

void Scores::postMergeStep() {
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
  totalNumberOfDecoys_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isDecoy)));
  totalNumberOfTargets_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isTarget)));
  targetDecoySizeRatio_ =
      totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  checkSeparationAndSetPi0();
}

double Scores::calcScore(const double* feat,
                         const std::vector<double>& w) const {
  std::size_t ix = FeatureNames::getNumFeatures();
  double score = w[ix];
  for (; ix--;) {
    score += feat[ix] * w[ix];
  }
  return score;
}

void Scores::scoreAndAddPSM(ScoreHolder& sh,
                            const std::vector<double>& rawWeights,
                            FeatureMemoryPool& featurePool) {
  const unsigned int numFeatures =
      static_cast<unsigned int>(FeatureNames::getNumFeatures());

  for (unsigned int j = 0; j < numFeatures; j++) {
    sh.score += sh.pPSM->features[j] * rawWeights[j];
  }
  sh.score += rawWeights[numFeatures];

  featurePool.deallocate(sh.pPSM->features);
  sh.pPSM->features = NULL;

  // DataSet::readPsm only returns TARGET or DECOY
  if (sh.isTarget()) {
    ++totalNumberOfTargets_;
  } else {
    ++totalNumberOfDecoys_;
  }
  scores_.push_back(sh);
}

void Scores::print(LabelType label, std::ostream& os) {
 std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  os << "PSMId\t";
  if (PSMDescription::hasSpectrumFileName()) {
    os << "filename\t";
  }
  os << "score\tq-value\tposterior_error_prob\tpeptide\tproteinIds";
  if (is_output_rt_) {
    os << "\trt\n";
  } else {
    os << "\n";
  }

  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->label == label) {
      PSMDescription* psm = scoreIt->pPSM;
      os << psm->getId() << '\t';
      const std::string& fileName = psm->getSpectrumFileName();
      if (!fileName.empty()) {
        os << fileName << '\t';
      }
      os << scoreIt->score << '\t' << scoreIt->q << '\t' << scoreIt->pep << '\t'
         << psm->getFullPeptide() << '\t';
      psm->printProteins(os);
      if (is_output_rt_) {
        os << '\t' << psm->getRetentionTime();
      }
      os << '\n';
    }
  }
}

void Scores::populateWithPSMs(SetHandler& setHandler) {
  scores_.clear();
  scores_.reserve(setHandler.getSizeFromLabel(LabelType::TARGET) +
                  setHandler.getSizeFromLabel(LabelType::DECOY));
  setHandler.populateScoresWithPSMs(scores_, LabelType::TARGET);
  setHandler.populateScoresWithPSMs(scores_, LabelType::DECOY);
  totalNumberOfTargets_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::TARGET));
  totalNumberOfDecoys_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::DECOY));
  targetDecoySizeRatio_ = (double)totalNumberOfTargets_ / totalNumberOfDecoys_;

  if (VERB > 1) {
    cerr << "Train/test set contains " << totalNumberOfTargets_
         << " positives and " << totalNumberOfDecoys_
         << " negatives, size ratio=" << targetDecoySizeRatio_
         << " and pi0=" << pi0_ << endl;
  }

  if (totalNumberOfTargets_ == 0) {
    ostringstream oss;
    oss << "Error: no target PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  if (totalNumberOfDecoys_ == 0) {
    ostringstream oss;
    oss << "Error: no decoy PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  // check for the minimum recommended number of positive and negative hits
  if (totalNumberOfTargets_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of positive samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
  if (totalNumberOfDecoys_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of negative samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
}

std::size_t Scores::getNumSpectra() const {
  std::vector<std::pair<unsigned int, unsigned int> > spectra;
  spectra.reserve(scores_.size());
  for (const ScoreHolder& sh : scores_) {
    spectra.emplace_back(sh.pPSM->specFileNr, sh.pPSM->scan);
  }
  std::sort(spectra.begin(), spectra.end());
  return static_cast<std::size_t>(
      std::unique(spectra.begin(), spectra.end()) - spectra.begin());
}

/**
 * Divides the PSMs from pin file into xval_fold cross-validation sets based on
 * their spectrum scan number. Each PSM is put in the test set of one fold and
 * in the training sets of the other xval_fold - 1 folds, and all PSMs of a
 * spectrum go to the same fold. With xval_fold = 1 the only training set and
 * the only test set both hold all PSMs.
 * @param train vector containing the training sets of PSMs
 * @param test vector containing the test sets of PSMs
 * @param xval_fold number of folds in train and test
 * @param decoyFractionTraining fraction of decoys to use for training, will
 * be 1 unless RESET algorithm is used
 */
void Scores::createXvalSetsBySpectrum(std::vector<Scores>& train,
                                      std::vector<Scores>& test,
                                      unsigned int xval_fold,
                                      double decoyFractionTraining,
                                      unsigned int decoysPerTarget) {
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold, Scores(usePi0_));
  test.resize(xval_fold, Scores(usePi0_));
  // The test folds are filled up to scores_.size() / xval_fold PSMs, rounded
  // up or down (remain, below), and a training set holds the PSMs of the
  // other test folds, or all PSMs with a single fold. Spectra with several
  // PSMs can make a fold a few PSMs larger.
  size_t approxTestSize = scores_.size() / xval_fold + 1;
  size_t approxTrainSize = (xval_fold > 1)
                               ? scores_.size() - scores_.size() / xval_fold
                               : scores_.size();
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].reserve(approxTrainSize);
    test[i].reserve(approxTestSize);
  }
  // remain keeps track of residual space available in each fold
  std::vector<int> remain(xval_fold);
  // set values for remain: initially each fold is assigned (tot number of
  // scores_ / tot number of folds)
  int fold = static_cast<int>(xval_fold), ix = static_cast<int>(scores_.size());
  while (fold--) {
    remain[static_cast<std::size_t>(fold)] = ix / (fold + 1);
    ix -= remain[static_cast<std::size_t>(fold)];
  }

  std::sort(scores_.begin(), scores_.end(), OrderScanHash());

  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  // put scores into the folds; choose a fold (at random) and change it only
  // when scores from a new spectra are encountered
  unsigned int previousSpectrum = scores_.begin()->pPSM->scan;
  size_t randIndex = PseudoRandom::lcg_rand() % xval_fold;
  for (std::vector<ScoreHolder>::iterator it = scores_.begin();
       it != scores_.end(); ++it) {
    const unsigned int curScan = (*it).pPSM->scan;
    ScoreHolder sh = (*it);

    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    if (previousSpectrum != curScan) {
      randIndex = PseudoRandom::lcg_rand() % xval_fold;
      // allow only indexes of folds that are non-full
      while (remain[randIndex] <= 0) {
        randIndex = PseudoRandom::lcg_rand() % xval_fold;
      }
    }

    // if we use multiple folds with RESET, assign 1-decoyFractionTraining as
    // pseudo targets.
    if (xval_fold > 1u && decoyFractionTraining < 1.0 &&
        sh.label == LabelType::DECOY &&
        PseudoRandom::lcg_uniform_rand() > decoyFractionTraining) {
      // From Algorithm S3 of the percolator-RESET supplementary material
      // decoyFractionTraining - the probability of assigning a decoy to the
      // training set
      sh.label = LabelType::PSEUDO_TARGET;
    }

    // insert
    for (unsigned int i = 0; i < xval_fold; ++i) {
      if (i == randIndex) {
        test[i].addScoreHolder(sh);
      }
      if (i != randIndex || xval_fold == 1) {
        train[i].addScoreHolder(sh);
      }
    }
    // update number of free position for used fold
    --remain[randIndex];
    // set previous spectrum to current one for next iteration
    previousSpectrum = curScan;
  }

  // without RESET: decoysPerTarget = 1 and decoyFractionTraining = 1, then
  // factor = 1 and trainNullTargetWinProb = 1/2 (and decoyFactor = 1)
  //
  // with RESET 1 decoy: decoysPerTarget = 1 and decoyFractionTraining = 0.5,
  // then factor = 1.5 and trainNullTargetWinProb = 3/4 (and decoyFactor = 3)
  //
  // with RESET 2 decoys: decoysPerTarget = 2 and decoyFractionTraining = 0.5,
  // then factor = 2 and trainNullTargetWinProb = 2/3 (and decoyFactor = 2)
  double factor =
      (decoysPerTarget + 1) - decoyFractionTraining * decoysPerTarget;
  double trainNullTargetWinProb = factor / (1.0 + decoysPerTarget);

  // calculate ratios of target over decoy for train and test set
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].recalculateSizes();
    test[i].recalculateSizes();

    train[i].setNullTargetWinProb(trainNullTargetWinProb);
    test[i].setNullTargetWinProb(trainNullTargetWinProb);
  }
}

void Scores::recalculateSizes() {
  totalNumberOfTargets_ = 0;
  totalNumberOfDecoys_ = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      ++totalNumberOfTargets_;
    } else {
      ++totalNumberOfDecoys_;
    }
  }
  targetDecoySizeRatio_ = totalNumberOfTargets_ / (double)totalNumberOfDecoys_;
}

/**
 * Linear rescaling of SVM scores such that q=fdr = 0 and the median decoy = -1.
 * @param fdr fdr-cutoff to use, typically --trainFDR
 * @param weights SVM weights used to compute SVM scores
 * @param singleBin whether these are the scores of the only cross-validation
 * bin; if they cannot be rescaled, the raw SVM scores and weights are then
 * kept with a warning instead of stopping, since a single bin needs no common
 * scale
 */
void Scores::normalizeScores(double fdr, std::vector<double>& weights,
                             bool singleBin) {
  unsigned int medianIndex = std::max(0u, totalNumberOfDecoys_ / 2u),
               decoys = 0u;

  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator it = scores_.begin();
  double fdrScore = it->score;
  double medianDecoyScore = fdrScore + 1.0;

  for (; it != scores_.end(); ++it) {
    if (it->q < fdr)
      fdrScore = it->score;
    if (it->isDecoy()) {
      if (++decoys == medianIndex) {
        medianDecoyScore = it->score;
        break;
      }
    }
  }

  double diff = fdrScore - medianDecoyScore;
  if (diff <= 0) {
    ostringstream oss;
    oss << "Error: median decoy score <= score at " << fdr * 100
        << "% FDR. Cannot rescale scores to merge cross validation bins, try "
        << "lowering --trainFDR.\n";
    if (NO_TERMINATE) {
      cerr << oss.str()
           << "No-terminate flag set: apply offset such that median "
           << "decoy has score -1, but skipping rescaling." << std::endl;
      diff = 1.0;
      fdrScore += 1.0;
    } else if (singleBin) {
      cerr << "Warning: median decoy score <= score at " << fdr * 100
           << "% FDR, the scores cannot be rescaled. Keeping the raw SVM "
           << "scores, as there is only one cross-validation bin to merge."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= fdrScore;
    scoreIt->score /= diff;
  }
  Normalizer::endScoreNormalizeWeights(weights, weights, fdrScore, diff);
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * and calculate q values for each PSM
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcScoresAndQvals(std::vector<double>& w,
                               double fdr,
                               bool skipDecoysPlusOne) {
  calcScores(w);
  return calcQvals(fdr, skipDecoysPlusOne);
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return 0 on successfull execution
 */
int Scores::calcScores(std::vector<double>& w) {
  std::size_t ix;
  const long numScores = static_cast<long>(scores_.size());
#pragma omp parallel for schedule(static) \
    if (numScores > 100000 && !omp_in_parallel())
  for (long i = 0; i < numScores; ++i) {
    scores_[i].score = calcScore(scores_[i].pPSM->features, w);
  }
  // already sorted when scored again with the same weights, as the test
  // sets are after training
  if (!isSortedByScore()) {
    sort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
  }
  if (VERB > 3) {
    if (scores_.size() >= 10) {
      cerr << "10 best scores and labels" << endl;
      for (ix = 0; ix < 10; ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
      cerr << "10 worst scores and labels" << endl;
      for (ix = scores_.size() - 10; ix < scores_.size(); ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
    } else {
      cerr << "Too few scores to display top and bottom PSMs ("
           << scores_.size() << " scores found)." << endl;
    }
  }
  return 0;
}

/**
 * Whether the PSMs are in strictly decreasing order of operator>, which
 * orders PSMs without NaN scores totally: then sorting them with
 * greater<ScoreHolder> would not change their order. A NaN score, or the
 * same PSM twice, gives false. Stops at the first pair out of order.
 */
bool Scores::isSortedByScore() const {
  for (std::size_t i = 1; i < scores_.size(); ++i) {
    if (!(scores_[i - 1] > scores_[i])) {
      return false;
    }
  }
  return true;
}

void Scores::getScoreLabelPairs(std::vector<pair<double, bool> >& combined) {
  combined.resize(scores_.size());
  transform(scores_.begin(), scores_.end(), combined.begin(),
            mem_fn(&ScoreHolder::toPair));
}

/**
 * Calculates the q-value for each psm in scores_: the q-value is the minimal
 * FDR of any set that includes the particular psm
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQvals(double fdr, bool skipDecoysPlusOne) {
  std::vector<pair<double, bool> > combined;
  return calcQvals(fdr, skipDecoysPlusOne, combined);
}

int Scores::calcQvals(double fdr, bool skipDecoysPlusOne,
                      std::vector<pair<double, bool> >& combined) {
  assert(totalNumberOfDecoys_ + totalNumberOfTargets_ == size());

  getScoreLabelPairs(combined);

  std::vector<double> qvals;
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne,
                                 nullTargetWinProb_);

  // set q-values and count number of positives
  // (countTargetsBelowFdr below has to give the same count)
  std::vector<double>::const_iterator qIt = qvals.begin();
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();

  int numPos = 0;
  for (; qIt != qvals.end(); ++qIt, ++scoreIt) {
    scoreIt->q = *qIt;
    if (scoreIt->q < fdr && scoreIt->isTarget())
      ++numPos;
  }

  return numPos;
}

/**
 * Counts the targets with q-value below fdr, as calcQvals does, for scores
 * listed in combined from best to worst. Within a group of tied scores all
 * PSMs get the same q-value, so the count does not depend on the order of
 * the PSMs within ties.
 * PosteriorEstimator::setNegative(true) must have been called.
 */
int Scores::countTargetsBelowFdr(
    const std::vector<pair<double, bool> >& combined,
    double fdr,
    bool skipDecoysPlusOne) const {
  std::vector<double> qvals;
  qvals.reserve(combined.size());
  PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne,
                                 nullTargetWinProb_);
  int numPos = 0;
  for (std::size_t ix = 0; ix < qvals.size(); ++ix) {
    if (qvals[ix] < fdr && combined[ix].second) {
      ++numPos;
    }
  }
  return numPos;
}

/**
 * Lists the feature rows and labels of the PSMs ordered by memory address,
 * so that scoring them streams through the feature storage instead of
 * following a pointer per PSM to scattered memory.
 */
void Scores::getFeatureRows(FeatureRows& rows) const {
  rows.resize(scores_.size());
  for (std::size_t ix = 0; ix < scores_.size(); ++ix) {
    assert((reinterpret_cast<std::uintptr_t>(scores_[ix].pPSM->features) &
            1u) == 0);
    rows[ix] = makeFeatureRow(scores_[ix].pPSM->features,
                              scores_[ix].isTarget());
  }
  // the PSMs have distinct feature rows, so the order is unique
  if (inParallelRegion()) {
    sequentialSort(rows.begin(), rows.end());
  } else {
    std::sort(rows.begin(), rows.end());
  }
}

namespace {

// The FDR estimate of PosteriorEstimator::getQValues without the mix-max
// correction (pi0 = 1), written as the same expression so that it gives
// the same bits.
inline double tdcFdrEstimate(int numDecoys, int numTargets,
                             double decoyFactor) {
  const double pi0 = 1.0, mixMaxCorrection = 0.0;
  return (numDecoys * pi0 + mixMaxCorrection) /
         (double)((std::max)(1, numTargets)) * decoyFactor;
}

// Sorts the PSMs ahead of a cut-off, best first, walks them as getQValues
// does and returns the number of targets with a q value below fdr. The q
// value of a PSM is the smallest estimate from its group of tied values
// on, so these are the targets up to the end of the last group whose
// estimate (at most 1) is below fdr.
template <class Better>
int countSortedTargetsBelowFdr(std::vector<pair<double, bool> >& psms,
                               Better better, double fdr,
                               int initialDecoyCount, double decoyFactor) {
  sequentialSort(psms.begin(), psms.end(),
                 [better](const pair<double, bool>& a,
                          const pair<double, bool>& b) {
                   return better(a.first, b.first);
                 });
  int numDecoys = initialDecoyCount, numTargets = 0, count = 0;
  const std::size_t n = psms.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (psms[i].second) {
      ++numTargets;
    } else {
      ++numDecoys;
    }
    if (i + 1 == n || psms[i].first != psms[i + 1].first) {
      double estimate = tdcFdrEstimate(numDecoys, numTargets, decoyFactor);
      if ((std::min)(estimate, 1.0) < fdr) {
        count = numTargets;
      }
    }
  }
  return count;
}

}  // namespace

bool Scores::countTargetsBySelection(const FeatureRows& rows,
                                     CountBuffers& buffers,
                                     double fdr,
                                     bool skipDecoysPlusOne,
                                     double nullTargetWinProb,
                                     int* highCount,
                                     int* lowCount) {
  const double decoyFactor = nullTargetWinProb / (1.0 - nullTargetWinProb);
  if (!(decoyFactor >= 0.0 &&
        decoyFactor <= std::numeric_limits<double>::max())) {
    return false;
  }
  const std::vector<double>& values = buffers.values;
  std::vector<double>& decoyValues = buffers.decoyValues;
  const std::size_t n = rows.size();
  assert(values.size() >= n);
  decoyValues.clear();
  int numTargets = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double value = values[i];
    if (value != value) {
      return false;  // NaN: a sort of all PSMs has no defined order
    }
    if (rowIsTarget(rows[i])) {
      ++numTargets;
    } else {
      decoyValues.push_back(value);
    }
  }
  const int numDecoys = static_cast<int>(decoyValues.size());
  const int initialDecoyCount = skipDecoysPlusOne ? 0 : 1;

  // The smallest number of decoys, cutOffDecoys, whose estimate with all
  // targets is not below fdr. Every cut-off with at least that many decoys
  // has at most all targets, and since IEEE division and multiplication are
  // monotone its estimate is not below fdr either. Such cut-offs cannot
  // give a q value below fdr, so only the PSMs with a better value than the
  // cutOffDecoys-th best decoy can count, and only they are sorted. If
  // there is no such number, all PSMs are sorted.
  int first = 0, last = numDecoys + 1;
  while (first < last) {
    const int mid = first + (last - first) / 2;
    double estimate =
        tdcFdrEstimate(initialDecoyCount + mid, numTargets, decoyFactor);
    if ((std::min)(estimate, 1.0) >= fdr) {
      last = mid;
    } else {
      first = mid + 1;
    }
  }
  const int cutOffDecoys = first;
  if (cutOffDecoys == 0) {  // no cut-off reaches q < fdr
    if (highCount) *highCount = 0;
    if (lowCount) *lowCount = 0;
    return true;
  }

  std::vector<pair<double, bool> >& high = buffers.high;
  std::vector<pair<double, bool> >& low = buffers.low;
  high.clear();
  low.clear();
  if (cutOffDecoys > numDecoys) {
    for (std::size_t i = 0; i < n; ++i) {
      const pair<double, bool> psm(values[i], rowIsTarget(rows[i]));
      if (highCount) high.push_back(psm);
      if (lowCount) low.push_back(psm);
    }
  } else {
    // the cutOffDecoys-th smallest and the cutOffDecoys-th largest decoy
    // value, at ascending ranks lowRank and highRank
    const int lowRank = cutOffDecoys - 1, highRank = numDecoys - cutOffDecoys;
    double lowCut = 0.0, highCut = 0.0;
    std::vector<double>::iterator dv = decoyValues.begin();
    if (lowCount) {
      sequentialNthElement(dv, dv + lowRank, decoyValues.end());
      lowCut = dv[lowRank];
    }
    if (highCount) {
      if (!lowCount) {
        sequentialNthElement(dv, dv + highRank, decoyValues.end());
      } else if (highRank > lowRank) {
        // the values from lowRank + 1 on are the larger ones
        sequentialNthElement(dv + lowRank + 1, dv + highRank,
                             decoyValues.end());
      } else if (highRank < lowRank) {
        sequentialNthElement(dv, dv + highRank, dv + lowRank);
      }
      highCut = dv[highRank];
    }
    for (std::size_t i = 0; i < n; ++i) {
      const double value = values[i];
      if (highCount && value > highCut) {
        high.push_back(pair<double, bool>(value, rowIsTarget(rows[i])));
      }
      if (lowCount && value < lowCut) {
        low.push_back(pair<double, bool>(value, rowIsTarget(rows[i])));
      }
    }
  }
  if (highCount) {
    *highCount = countSortedTargetsBelowFdr(high, std::greater<double>(), fdr,
                                            initialDecoyCount, decoyFactor);
  }
  if (lowCount) {
    *lowCount = countSortedTargetsBelowFdr(low, std::less<double>(), fdr,
                                           initialDecoyCount, decoyFactor);
  }
  return true;
}

/**
 * Returns the number of targets that calcScoresAndQvals(w, fdr,
 * skipDecoysPlusOne) would report, without changing this object, so that it
 * can be used concurrently. Only scores and labels are sorted, which is much
 * cheaper than sorting the ScoreHolders, and without the mix-max correction
 * (pi0 = 1) only those ahead of the cut-off decoy of countTargetsBySelection.
 * @param rows the PSMs of this set, from getFeatureRows
 * @param buffers the buffers of the calling thread
 */
int Scores::countTargetsBelowFdr(const FeatureRows& rows,
                                 const std::vector<double>& w,
                                 double fdr,
                                 bool skipDecoysPlusOne,
                                 CountBuffers& buffers) const {
  // same summation order as calcScore, four rows at a time so that their
  // additions overlap
  const std::size_t numFeatures = FeatureNames::getNumFeatures();
  const std::size_t numRows = rows.size();
  std::vector<double>& scores = buffers.values;
  scores.resize(numRows);
  std::size_t i = 0;
  for (; i + 4 <= numRows; i += 4) {
    const double* f0 = rowFeatures(rows[i]);
    const double* f1 = rowFeatures(rows[i + 1]);
    const double* f2 = rowFeatures(rows[i + 2]);
    const double* f3 = rowFeatures(rows[i + 3]);
    double s0 = w[numFeatures], s1 = w[numFeatures], s2 = w[numFeatures],
           s3 = w[numFeatures];
    for (std::size_t ix = numFeatures; ix--;) {
      s0 += f0[ix] * w[ix];
      s1 += f1[ix] * w[ix];
      s2 += f2[ix] * w[ix];
      s3 += f3[ix] * w[ix];
    }
    scores[i] = s0;
    scores[i + 1] = s1;
    scores[i + 2] = s2;
    scores[i + 3] = s3;
  }
  for (; i < numRows; ++i) {
    scores[i] = calcScore(rowFeatures(rows[i]), w);
  }
  int count = 0;
  if (pi0_ == 1.0 &&
      countTargetsBySelection(rows, buffers, fdr, skipDecoysPlusOne,
                              nullTargetWinProb_, &count, nullptr)) {
    return count;
  }
  // mix-max counts and NaN scores: sort all PSMs
  std::vector<pair<double, bool> > combined(numRows);
  for (i = 0; i < numRows; ++i) {
    combined[i] = pair<double, bool>(scores[i], rowIsTarget(rows[i]));
  }
  std::sort(combined.begin(), combined.end(),
            [](const pair<double, bool>& a, const pair<double, bool>& b) {
              return a.first > b.first;
            });
  return countTargetsBelowFdr(combined, fdr, skipDecoysPlusOne);
}

/**
 * Starts a new training set in data with the decoys, in the order of this
 * set. generatePositiveTrainingSet then appends the positives.
 */
void Scores::generateNegativeTrainingSet(AlgIn& data, const double cneg) {
  data.clear();
  // room for the largest possible training set; only touched when used
  data.rows.reserve(scores_.size() * static_cast<size_t>(data.n));
  data.Y.reserve(scores_.size());
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isDecoy()) {
      data.add(scoreIt->pPSM->features, -1);
    }
  }
  data.negatives = data.m;
}

/**
 * Appends the targets with q-value at most fdr to the training set, in the
 * order of this set, or with trainBestPositive only the best-scoring one of
 * each spectrum, by decreasing score. Does not change this set.
 */
void Scores::generatePositiveTrainingSet(AlgIn& data,
                                         const double fdr,
                                         const double cpos,
                                         const bool trainBestPositive) const {
  int p = 0;
  if (trainBestPositive) {
    std::vector<ScoreHolder> best;
    for (const ScoreHolder& sh : scores_) {
      if (sh.isTarget() && sh.q <= fdr) best.push_back(sh);
    }
    std::sort(best.begin(), best.end(), OrderScanLabel());
    best.erase(std::unique(best.begin(), best.end(), UniqueScanLabel()),
               best.end());
    std::sort(best.begin(), best.end(), greater<ScoreHolder>());
    for (const ScoreHolder& sh : best) {
      data.add(sh.pPSM->features, 1);
      ++p;
    }
  } else {
    for (const ScoreHolder& sh : scores_) {
      if (sh.isTarget() && sh.q <= fdr) {
        data.add(sh.pPSM->features, 1);
        ++p;
      }
    }
  }
  data.positives = p;
}

void Scores::weedOutRedundant() {
  std::map<std::string, unsigned int> peptideSpecCounts;
  double specCountQvalThreshold = -1.0;
  weedOutRedundant(peptideSpecCounts, specCountQvalThreshold);
}

/**
 * Routine that sees to that only unique peptides are kept (used for analysis
 * on peptide-fdr rather than psm-fdr)
 */
void Scores::weedOutRedundant(
    std::map<std::string, unsigned int>& peptideSpecCounts,
    double specCountQvalThreshold) {
  // lexicographically order the scores_ (based on peptides names,labels and
  // scores)
  std::sort(scores_.begin(), scores_.end(), lexicOrderProb());

  // the views point into the PSMs, which stay in place
  std::string_view previousPeptide = "";
  LabelType previousLabel = LabelType::UNDEFINED;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx) {
    std::string_view currentPeptide =
        scores_[idx].pPSM->getPeptideSequenceView();
    LabelType currentLabel = scores_[idx].label;
    if (currentPeptide != previousPeptide || currentLabel != previousLabel) {
      // insert as a new score
      scores_[lastWrittenIdx++] = scores_[idx];
      previousPeptide = currentPeptide;
      previousLabel = currentLabel;
    }
    if (specCountQvalThreshold > 0.0 &&
        scores_[idx].q < specCountQvalThreshold) {
      ++peptideSpecCounts[std::string(currentPeptide)];
    }
  }
  scores_.resize(lastWrittenIdx);
  postMergeStep();
}

/**
 * Routine that sees to that only unique spectra are kept for TDC
 */
void Scores::weedOutRedundantTDC() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());
  scores_.erase(
      std::unique(scores_.begin(), scores_.end(), UniqueScanMassCharge()),
      scores_.end());

  /* does not actually release memory because of memory fragmentation
  double previousExpMass = 0.0;
  unsigned int previousScan = 0u;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx){
    double currentExpMass = scores_.at(idx).pPSM->expMass;
    int currentScan = scores_.at(idx).pPSM->scan;
    if (currentExpMass != previousExpMass || currentScan != previousScan) {
      // insert as a new score
      scores_.at(lastWrittenIdx++).swap(scores_.at(idx));
      previousScan = currentScan;
      previousExpMass = currentExpMass;
    } else {
      PSMDescription::deletePtr(scores_.at(idx).pPSM);
    }
  }
  scores_.resize(lastWrittenIdx);
  */
  postMergeStep();
}

int Scores::getInitDirection(const double initialSelectionFdr,
                             std::vector<double>& direction) {
  int bestPositives = -1;
  int bestFeature = -1;
  bool lowBest = false;

  // for determining the initial direction, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true;

  // Number of positives when ordering by each feature, with low values
  // first (i = 0) or high values first (i = 1). They only depend on the
  // feature values and labels, so the features are evaluated in parallel.
  // Each column is read along the feature rows in memory order, and
  // without the mix-max correction both counts come from one selection
  // (countTargetsBySelection).
  const int numFeatures = static_cast<int>(FeatureNames::getNumFeatures());
  std::vector<int> positivesPerDirection(2 * static_cast<size_t>(numFeatures));
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  FeatureRows rows;
  getFeatureRows(rows);
  const std::size_t numRows = rows.size();
#pragma omp parallel
  {
    CountBuffers buffers;  // one per thread, for all its features
    buffers.values.resize(numRows);
#pragma omp for schedule(dynamic, 1)
    for (int featNo = 0; featNo < numFeatures; featNo++) {
      for (std::size_t ix = 0; ix < numRows; ++ix) {
        buffers.values[ix] = rowFeatures(rows[ix])[featNo];
      }
      int lowCount = 0, highCount = 0;
      if (pi0_ == 1.0 &&
          countTargetsBySelection(rows, buffers, initialSelectionFdr,
                                  skipDecoysPlusOne, nullTargetWinProb_,
                                  &highCount, &lowCount)) {
        positivesPerDirection[2 * featNo] = lowCount;
        positivesPerDirection[2 * featNo + 1] = highCount;
        continue;
      }
      // mix-max counts and NaN values: sort all PSMs
      std::vector<pair<double, bool> > combined(scores_.size());
      for (std::size_t ix = 0; ix < scores_.size(); ++ix) {
        combined[ix] = pair<double, bool>(scores_[ix].pPSM->features[featNo],
                                          scores_[ix].isTarget());
      }
      std::sort(combined.begin(), combined.end(),
                [](const pair<double, bool>& a, const pair<double, bool>& b) {
                  return a.first < b.first;
                });
      positivesPerDirection[2 * featNo] = countTargetsBelowFdr(
          combined, initialSelectionFdr, skipDecoysPlusOne);
      std::reverse(combined.begin(), combined.end());
      positivesPerDirection[2 * featNo + 1] = countTargetsBelowFdr(
          combined, initialSelectionFdr, skipDecoysPlusOne);
    }
  }

  for (int featNo = 0; featNo < numFeatures; featNo++) {
    for (int i = 0; i < 2; i++) {
      int positives = positivesPerDirection[2 * featNo + i];
      if (positives > bestPositives) {
        bestPositives = positives;
        bestFeature = featNo;
        lowBest = (i == 0);
      }
    }
  }
  for (std::size_t ix = FeatureNames::getNumFeatures(); ix--;) {
    direction[ix] = 0;
  }

  if (bestPositives <= 0) {
    ostringstream oss;
    oss << "Error in the input data: cannot find an initial direction with "
        << "positive training examples. "
        << "Consider setting/raising the initial training FDR threshold "
           "(--train-fdr-initial)."
        << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str();
      std::cerr << "No-terminate flag set: setting initial direction to the "
                << "first feature and ignoring the error." << std::endl;
      bestFeature = 0;
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }

  if (bestFeature >= 0) {
    direction[static_cast<std::size_t>(bestFeature)] = (lowBest ? -1 : 1);
  }

  if (VERB > 1) {
    cerr << "Selected feature " << bestFeature + 1
         << " as initial direction. Could separate " << bestPositives
         << " training set positives with q<" << initialSelectionFdr
         << " in that direction." << endl;
  }
  return bestPositives;
}

void Scores::checkSeparationAndSetPi0() {
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);

  std::vector<double> pvals;
  PosteriorEstimator::getPValues(combined, pvals);

  pi0_ = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvals);
  if (tooGoodSeparation) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (usePi0_) {
        std::cerr
            << "No-terminate flag set: setting pi0 = 1 and ignoring error."
            << std::endl;
      } else {
        std::cerr << "No-terminate flag set: ignoring error." << std::endl;
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  } else if (usePi0_) {
    pi0_ = PosteriorEstimator::estimatePi0(pvals);
  }
}

void Scores::calcPep(const bool spline, const bool interp, const bool pava,
                     const bool tdc) {
    if (!spline && tdc) {
        // Fit the decoy rate p as a monotone function of the score and
        // report the odds p / (1 - p) as PEP, for targets and decoys alike.
        std::vector<double> is_decoy, sc;
        is_decoy.reserve(scores_.size());
        sc.reserve(scores_.size());
        for (const auto& sh : scores_) {
            is_decoy.push_back(sh.isTarget() ? 0.0 : 1.0);
            sc.push_back(sh.score);
        }
        // With separate target and decoy searches (mix-max) the decoys
        // sample the null of all targets, so the odds are scaled by pi0 and
        // by the ratio of targets to decoys.
        double oddsFactor = 1.0;
        if (usePi0_) {
            double numTargets = 0.0, numDecoys = 0.0;
            for (double d : is_decoy) (d > 0.5 ? numDecoys : numTargets) += 1.0;
            oddsFactor = pi0_ * numTargets / std::max(1.0, numDecoys);
        }
        InferPEP reg(!pava);
        auto peps = interp ? reg.tdc_to_pep(is_decoy, sc, oddsFactor)
                           : reg.tdc_to_pep(is_decoy, {}, oddsFactor);
        for (size_t ix = 0; ix < scores_.size(); ++ix) {
            scores_[ix].pep = peps[ix];
        }
    } else if (!spline) {
        // Derive PEPs from q-values so that average PEP ≈ q at every
        // threshold.  Use PAVA when --pava-pep is set, otherwise I-spline
        // for smooth output.
        const bool use_ispline = !pava;
        std::vector<double> target_q, sc;
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                target_q.push_back(sh.q);
                sc.push_back(sh.score);
            }
        }
        InferPEP reg(use_ispline);
        auto target_pep = interp
                            ? reg.qns_to_pep(target_q, sc)
                            : reg.q_to_pep(target_q);
        // Move PEPs to scoreholders. The PEPs are only defined for target,
        // We use interpolation for decoys.
        // Add elements avoiding overflow problems if last sh is a decoy
        target_pep.push_back(1.0); target_q.push_back(1.0);
        auto it_pep = target_pep.begin();
        auto it_q = target_q.begin();
        double l_q(0.0), l_pep(0.0);
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                sh.pep = *it_pep;
                // remember last (l_) pep and q for interpolation
                l_pep = *it_pep;
                l_q = *it_q;
                it_pep++; it_q++;
            } else {
                double pep = reg.interpolate(sh.q,l_q,*it_q,l_pep,*it_pep);
                sh.pep = std::max(0.0, std::min(1.0, pep));
            }
        }
    } else {
        std::vector<pair<double, bool> > combined;
        getScoreLabelPairs(combined);

        std::vector<double> peps;
        // Logistic regression on the data
        PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, true);
        for (size_t ix = 0; ix < scores_.size(); ix++) {
            scores_[ix].pep = peps[ix];
        }
    }
}

unsigned Scores::getQvaluesBelowLevel(double level) {
  unsigned hits = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() && scoreIt->q < level) {
      hits++;
    }
  }
  return hits;
}

void Scores::setUsePi0(bool usePi0) {
  usePi0_ = usePi0;
}
