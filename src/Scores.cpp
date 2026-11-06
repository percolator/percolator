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
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DataSet.h"
#include "Globals.h"
#include "MassHandler.h"
#include "Normalizer.h"
#include "ParallelSort.h"
#include "PosteriorEstimator.h"
#include "Scores.h"
#include "SetHandler.h"
#include "SequentialSort.h"
#include "ssl.h"
#include "IsotonicPEP.h"
#include "NumberFormatter.h"

void Scores::merge(std::vector<Scores>& sv,
                   double fdr,
                   bool skipNormalizeScores,
                   std::vector<std::vector<double> >& all_w,
                   MergeNormalization scheme,
                   const std::vector<ScoreScale>& binScales) {
  const bool fromTrainingSets = scheme == MergeNormalization::TRAIN ||
                                scheme == MergeNormalization::TRAIN_DECOYS;
  if (!skipNormalizeScores && fromTrainingSets &&
      binScales.size() != sv.size()) {
    throw MyException(
        "Error: the training-set rescaling of the cross-validation bins is "
        "missing.\n");
  }
  scores_.clear();
  size_t numScores = 0;
  for (const Scores& binScores : sv) {
    numScores += binScores.scores_.size();
  }
  scores_.reserve(numScores);
  std::vector<Scores>::iterator cvBinScores = sv.begin();
  std::vector<std::vector<double> >::iterator weights = all_w.begin();
  for (std::size_t bin = 0; cvBinScores != sv.end();
       cvBinScores++, weights++, ++bin) {
    // the test sets are usually sorted already, by calcScores
    if (!cvBinScores->isSortedByScore()) {
      parallelSort(cvBinScores->begin(), cvBinScores->end(),
                   greater<ScoreHolder>());
    }
    // the q values of a bin are only needed for its own anchors; the merged
    // list gets new ones
    if (scheme == MergeNormalization::TEST) {
      cvBinScores->checkSeparationAndSetPi0();
      cvBinScores->calcQvals(fdr);
    }
    if (!skipNormalizeScores) {
      if (scheme == MergeNormalization::TEST) {
        cvBinScores->normalizeScores(fdr, *weights, sv.size() == 1u);
      } else if (fromTrainingSets) {
        cvBinScores->applyScoreScale(binScales[bin], *weights,
                                     sv.size() == 1u, scheme, fdr);
      }
    }
    copy(cvBinScores->begin(), cvBinScores->end(), back_inserter(scores_));
    // release the bin before the sort of all PSMs
    cvBinScores->reset();
  }
  postMergeStep();
}

void Scores::postMergeStep() {
  parallelSort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
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

  const NumberFormatter format(os);
  if (format.sameAsStream()) {
    // the lines are built in a buffer, which is written in large blocks
    std::string buffer;
    for (; scoreIt != scores_.end(); ++scoreIt) {
      if (scoreIt->label != label) continue;
      PSMDescription* psm = scoreIt->pPSM;
      buffer.append(psm->getId());
      buffer += '\t';
      const std::string& fileName = psm->getSpectrumFileName();
      if (!fileName.empty()) {
        buffer.append(fileName);
        buffer += '\t';
      }
      format.append(buffer, scoreIt->score);
      buffer += '\t';
      format.append(buffer, scoreIt->q);
      buffer += '\t';
      format.append(buffer, scoreIt->pep);
      buffer += '\t';
      buffer.append(psm->getFullPeptide());
      buffer += '\t';
      psm->appendProteins(buffer);
      if (is_output_rt_) {
        buffer += '\t';
        format.append(buffer, psm->getRetentionTime());
      }
      buffer += '\n';
      if (buffer.size() >= NumberFormatter::kWriteBlockSize) {
        os.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        buffer.clear();
      }
    }
    os.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    return;
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
 * the only test set both hold all PSMs. With at least as many spectra as
 * folds, every test set holds at least one spectrum; a spectrum ends where
 * the scan number changes in the order of the split.
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

  parallelSort(scores_.begin(), scores_.end(), OrderScanHash());

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

  // The folds are filled with whole spectra while remain counts PSMs, so a
  // fold can take more PSMs than its share, and with uneven spectra or many
  // folds the spectra could run out before every fold has one. Once the
  // spectra left are as many as the folds that are still empty, each of them
  // goes to an empty fold. A draw is only changed where it would have left a
  // fold empty: a spectrum that goes to a non-empty fold at that point
  // leaves more empty folds than spectra. The spectra are counted as the
  // loop below delimits them, by a change of the scan number between
  // neighbours.
  std::size_t spectraLeft = 1u;
  for (std::size_t i = 1u; i < scores_.size(); ++i) {
    if (scores_[i].pPSM->scan != scores_[i - 1u].pPSM->scan) {
      ++spectraLeft;
    }
  }
  std::vector<bool> foldIsEmpty(xval_fold, true);
  std::size_t numEmptyFolds = xval_fold;

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
    if (it == scores_.begin() || previousSpectrum != curScan) {
      // spectraLeft counts this spectrum and the ones after it
      if (spectraLeft == numEmptyFolds && !foldIsEmpty[randIndex]) {
        // the first empty fold, without a random draw
        randIndex = 0u;
        while (!foldIsEmpty[randIndex]) {
          ++randIndex;
        }
      }
      if (foldIsEmpty[randIndex]) {
        foldIsEmpty[randIndex] = false;
        --numEmptyFolds;
      }
      --spectraLeft;
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
  applyScoreScale(getScoreScale(MergeNormalization::TEST, fdr), weights,
                  singleBin, MergeNormalization::TEST, fdr);
}

ScoreScale Scores::getScoreScale(MergeNormalization scheme,
                                 double fdr) const {
  ScoreScale s;
  if (scheme == MergeNormalization::NONE) {
    s.scale = 1.0;
    return s;
  }
  if (scheme == MergeNormalization::TRAIN_DECOYS) {
    // mean and sample standard deviation of the decoy scores, in two passes
    double sum = 0.0;
    std::size_t numDecoys = 0u;
    for (const ScoreHolder& sh : scores_) {
      if (sh.isDecoy()) {
        sum += sh.score;
        ++numDecoys;
      }
    }
    if (numDecoys == 0u) {
      return s;
    }
    const double mean = sum / static_cast<double>(numDecoys);
    double sumSquares = 0.0;
    for (const ScoreHolder& sh : scores_) {
      if (sh.isDecoy()) {
        sumSquares += (sh.score - mean) * (sh.score - mean);
      }
    }
    s.offset = mean;
    s.scale = numDecoys > 1u ? std::sqrt(sumSquares /
                                         static_cast<double>(numDecoys - 1u))
                             : 0.0;
    return s;
  }
  // TEST and TRAIN: the score of the last PSM with q < fdr above the median
  // decoy, or the top score when there is none, and the median decoy; with
  // fewer than two decoys there is no median decoy and the scale is negative
  if (scores_.empty()) {
    return s;
  }
  unsigned int medianIndex = std::max(0u, totalNumberOfDecoys_ / 2u),
               decoys = 0u;
  std::vector<ScoreHolder>::const_iterator it = scores_.begin();
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
  s.offset = fdrScore;
  s.scale = fdrScore - medianDecoyScore;
  return s;
}

ScoreScale Scores::calcScoreScale(std::vector<double>& w,
                                  MergeNormalization scheme,
                                  double fdr) {
  calcScores(w);
  if (scheme == MergeNormalization::TEST ||
      scheme == MergeNormalization::TRAIN) {
    checkSeparationAndSetPi0(false);
    calcQvals(fdr);
  }
  return getScoreScale(scheme, fdr);
}

void Scores::applyScoreScale(const ScoreScale& s,
                             std::vector<double>& weights,
                             bool singleBin,
                             MergeNormalization scheme,
                             double fdr) {
  double offset = s.offset;
  double scale = s.scale;
  if (scale <= 0) {
    const bool decoySpread = scheme == MergeNormalization::TRAIN_DECOYS;
    auto describe = [&](std::ostream& os) {
      if (decoySpread) {
        os << "the decoys of a training set have no spread in their scores";
      } else {
        os << "median decoy score <= score at " << fdr * 100 << "% FDR";
        if (scheme == MergeNormalization::TRAIN) {
          os << " in a training set";
        }
      }
    };
    ostringstream oss;
    oss << "Error: ";
    describe(oss);
    oss << ". Cannot rescale scores to merge cross validation bins, try "
        << (decoySpread ? "fewer folds (--num-folds)" : "lowering --trainFDR")
        << ".\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (decoySpread) {
        cerr << "No-terminate flag set: subtract the mean decoy score, but "
             << "skipping rescaling." << std::endl;
      } else {
        cerr << "No-terminate flag set: apply offset such that median "
             << "decoy has score -1, but skipping rescaling." << std::endl;
        offset += 1.0;
      }
      scale = 1.0;
    } else if (singleBin) {
      cerr << "Warning: ";
      describe(cerr);
      cerr << ", the scores cannot be rescaled. Keeping the raw SVM "
           << "scores, as there is only one cross-validation bin to merge."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= offset;
    scoreIt->score /= scale;
  }
  Normalizer::endScoreNormalizeWeights(weights, weights, offset, scale);
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
    parallelSort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
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
 * FDR of any set that includes the particular psm. The PSMs must be sorted
 * by decreasing score.
 *
 * Gives the same q values as PosteriorEstimator::getQValues on the (score,
 * isTarget) pairs of scores_ with setNegative(true), but writes them to the
 * PSMs directly instead of building the pairs, the q values and the mix-max
 * counts in arrays of their own.
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQvals(double fdr, bool skipDecoysPlusOne) {
  assert(totalNumberOfDecoys_ + totalNumberOfTargets_ == size());

  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  const std::size_t numScores = scores_.size();
  if (numScores == 0) {
    return 0;
  }
  const double pi0 = pi0_;
  const double decoyFactor = nullTargetWinProb_ / (1.0 - nullTargetWinProb_);

  // The mix-max correction of a tie group with decoys uses N_{w<=z} and
  // N_{z<=z}, the numbers of targets and decoys that score at most as much as
  // the (n_z_ge_w - 1)-th best decoy, which getQValues reads from the arrays
  // of getMixMaxCounts. They follow from the totals and the counts above
  // the group.
  int numTargets = 0, numDecoys = 0;
  if (pi0 < 1.0) {
    for (const ScoreHolder& sh : scores_) {
      (sh.isTarget() ? numTargets : numDecoys) += 1;
    }
  }
  // the counts at the last tie group with decoys
  int lastCntW = 0, lastCntZ = 0;

  double estPx_lt_zj = 0.0;
  double E_f1_mod_run_tot = 0.0;
  double fdrValue = 0.0;

  int n_z_ge_w = 1, n_w_ge_w = 0;  // N_{z>=w} and N_{w>=w}
  if (skipDecoysPlusOne) n_z_ge_w = 0;
  int numDecoysAbove = 0;  // decoys ranked at or above the current PSM

  int decoyQueue = 0, targetQueue = 0;  // handles ties
  std::size_t groupBegin = 0;
  for (std::size_t ix = 0; ix < numScores; ++ix) {
    if (scores_[ix].isTarget()) {
      ++n_w_ge_w;  // target PSM
      ++targetQueue;
    } else {
      ++n_z_ge_w;  // decoy PSM
      ++numDecoysAbove;
      ++decoyQueue;
    }

    // handles ties
    if (ix + 1 == numScores || scores_[ix].score != scores_[ix + 1].score) {
      if (pi0 < 1.0 && decoyQueue > 0) {
        const int cntW = numTargets - (n_w_ge_w - targetQueue);
        const int cntZ = numDecoys - (numDecoysAbove - decoyQueue);
        // rank of the decoy whose counts are used, 1 for the best one
        const int rank = n_z_ge_w - 1;
        int cnt_w = cntW, cnt_z = cntZ;
        if (rank <= numDecoysAbove - decoyQueue) {
          // with skipDecoysPlusOne the decoy ranked just above a group with
          // a single decoy, i.e. in the last group with decoys
          if (rank < 1) {
            throw std::out_of_range(
                "Scores::calcQvals: no mix-max counts above the best decoy");
          }
          cnt_w = lastCntW;
          cnt_z = lastCntZ;
        }
        lastCntW = cntW;
        lastCntZ = cntZ;
        estPx_lt_zj = (double)(cnt_w - pi0 * cnt_z) / ((1.0 - pi0) * cnt_z);
        estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
        estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
        E_f1_mod_run_tot += decoyQueue * estPx_lt_zj * (1.0 - pi0);
        if (VERB > 4) {
          std::cerr << "Mix-max num negatives correction: "
                    << (1.0 - pi0) * n_z_ge_w << " vs. " << E_f1_mod_run_tot
                    << std::endl;
        }
      }

      fdrValue = (n_z_ge_w * pi0 + E_f1_mod_run_tot) /
                 (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
      const double groupFdr = (std::min)(fdrValue, 1.0);
      for (std::size_t i = groupBegin; i <= ix; ++i) {
        scores_[i].q = groupFdr;
      }
      groupBegin = ix + 1;
      decoyQueue = 0;
      targetQueue = 0;
    }
  }

  // Convert the FDRs into q-values, from the bottom up, and count the
  // targets below fdr (countTargetsBelowFdr has to give the same count).
  int numPos = 0;
  double minFdr = scores_[numScores - 1].q;
  for (std::size_t ix = numScores; ix--;) {
    ScoreHolder& sh = scores_[ix];
    // the same comparison as mymin in PosteriorEstimator::getQValues
    minFdr = minFdr > sh.q ? sh.q : minFdr;
    sh.q = minFdr;
    if (sh.q < fdr && sh.isTarget()) ++numPos;
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
  // the PSMs have distinct feature rows, so the order is unique;
  // parallelSort sorts on the calling thread inside a parallel region
  parallelSort(rows.begin(), rows.end(), std::less<std::uintptr_t>());
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
void Scores::generateNegativeTrainingSet(AlgIn& data, const double /*cneg*/) {
  data.clear();
  // room for the largest possible training set; only touched when used
  data.reserve(scores_.size());
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
                                         const double /*cpos*/,
                                         const bool trainBestPositive) const {
  int p = 0;
  if (trainBestPositive) {
    std::vector<ScoreHolder> best;
    for (const ScoreHolder& sh : scores_) {
      if (sh.isTarget() && sh.q <= fdr) best.push_back(sh);
    }
    // with nested cross validation called in a parallel loop, where
    // parallelSort sorts on the calling thread
    parallelSort(best.begin(), best.end(), OrderScanLabel());
    best.erase(std::unique(best.begin(), best.end(), UniqueScanLabel()),
               best.end());
    parallelSort(best.begin(), best.end(), greater<ScoreHolder>());
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
  parallelSort(scores_.begin(), scores_.end(), lexicOrderProb());

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
  parallelSort(scores_.begin(), scores_.end(), OrderScanMassCharge());
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
  // the bias too: after a direction reset, direction holds the trained
  // weights of the fold
  std::fill(direction.begin(), direction.end(), 0.0);

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

/**
 * The p values of the targets, in ascending order, from the decoys that
 * score at least as much: the same as PosteriorEstimator::getPValues on the
 * (score, isTarget) pairs of scores_, without building the pairs. The PSMs
 * must be sorted by decreasing score.
 */
void Scores::getTargetPValues(std::vector<double>& pvals) const {
  std::size_t numTargets = 0;
  for (const ScoreHolder& sh : scores_) {
    if (sh.isTarget()) ++numTargets;
  }
  pvals.clear();
  pvals.reserve(numTargets);
  std::size_t nDecoys = 1, posSame = 0, negSame = 0;
  const std::size_t numScores = scores_.size();
  for (std::size_t ix = 0; ix < numScores; ++ix) {
    if (scores_[ix].isTarget()) {
      ++posSame;
    } else {
      ++negSame;
    }
    if (ix + 1 == numScores || scores_[ix].score != scores_[ix + 1].score) {
      for (std::size_t i = 0; i < posSame; ++i) {
        pvals.push_back(static_cast<double>(nDecoys) +
                        static_cast<double>(negSame * (i + 1)) /
                            (double)(posSame + 1));
      }
      nDecoys += negSame;
      negSame = 0;
      posSame = 0;
    }
  }
  for (double& p : pvals) {
    p = p / static_cast<double>(nDecoys);
  }
}

void Scores::checkSeparationAndSetPi0(bool tooGoodSeparationIsError) {
  std::vector<double> pvals;
  getTargetPValues(pvals);

  pi0_ = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvals);
  if (tooGoodSeparation && !tooGoodSeparationIsError) {
    return;
  }
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
    // Expected number of incorrect targets per decoy, as in the q values: 1
    // in target-decoy competition, 2 with RESET, which keeps only the half
    // of the decoys that were not used for training.
    const double decoyFactor = nullTargetWinProb_ / (1.0 - nullTargetWinProb_);
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
        oddsFactor *= decoyFactor;
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
        PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, true,
                                        decoyFactor);
        for (size_t ix = 0; ix < scores_.size(); ix++) {
            scores_[ix].pep = peps[ix];
        }
    }
}

void Scores::setUsePi0(bool usePi0) {
  usePi0_ = usePi0;
}
