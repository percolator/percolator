/*
 * Copyright 2021 Brian Raiter <breadbox@muppetlabs.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

/*
 * Unit tests for the Scores and ScoreHolder classes.
 */

#include <gtest/gtest.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "SetHandler.h"
#include "DataSet.h"
#include "PosteriorEstimator.h"
#include "Scores.h"

namespace {

// Runs OpenMP code, such as GCC's parallel-mode sorts and partial sums in
// the reference counts, on one thread while it exists: thousands of small
// parallel calls take minutes on a loaded machine.
class OneThread {
  public:
    OneThread() {
#ifdef _OPENMP
        numThreads_ = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
    }
    ~OneThread() {
#ifdef _OPENMP
        omp_set_num_threads(numThreads_);
#endif
    }
  private:
    int numThreads_ = 1;
};

}  // namespace

// Some strings in alphabetical order.
static std::string const psmNames[] = { "ABC", "DEF", "GHI", "JKL", "MNO" };

class ScoreHolderTest : public ::testing::Test {
  protected:
    bool checkOrder(std::vector<ScoreHolder> const *scores, ...);
    // a PSM with the given peptide, owned by the fixture
    PSMDescription* newPsm(const std::string& peptide) {
        PSMDescription* psm = psmStore.add(PSMDescription());
        psmStore.storeStrings(*psm, "", peptide, {});
        return psm;
    }
    PsmStore psmStore;
};

// A simple function that verifies that the scores in a vector are in
// a specific order. (The function is variadic so that the caller
// doesn't need to declare an extra array.)
bool ScoreHolderTest::checkOrder(std::vector<ScoreHolder> const *scores, ...)
{
    va_list values;
    va_start(values, scores);
    for (std::vector<ScoreHolder>::const_iterator it = scores->begin() ;
             it != scores->end() ; ++it) {
        double score = va_arg(values, double);
        if (it->score != score)
            return false;
    }
    va_end(values);
    return true;
}

// Verify the ScoreHolder's ordering functions.
TEST_F(ScoreHolderTest, CheckOrderingFunctions)
{
    std::vector<ScoreHolder> scores;

    // peptides are assigned in alphabetic order
    // scan values are assigned [ 2, 3, 4, 0, 1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = newPsm(psmNames[i]);
        pPSM->scan = (i + 2) % 5;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), lexicOrderProb());
    ASSERT_TRUE(checkOrder(&scores, 1.0, 2.0, 3.0, 4.0, 5.0));
    std::sort(scores.begin(), scores.end(), OrderScanMassCharge());
    ASSERT_TRUE(checkOrder(&scores, 4.0, 5.0, 1.0, 2.0, 3.0));
    
    // specFileNr values are assigned [ 3, 3, 4, 4, 5]
    // scan values are assigned [ 0, 0, 1, 1, 2]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = newPsm(psmNames[i]);
        pPSM->specFileNr = 3 + i / 2;
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanHash());
    for (int i = 0 ; i < 5 ; ++i) {
        std::cerr << scores.at(i).score << std::endl;
    }
    ASSERT_TRUE(checkOrder(&scores, 3.0, 4.0, 1.0, 2.0, 5.0));

    // specFileNr values are assigned [ 0, 0, 1, 1, 2]
    // scan values are assigned [ 2, 3, 4, 0, 1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = newPsm(psmNames[i]);
        pPSM->specFileNr = i / 2;
        pPSM->scan = (i + 2) % 5;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanMassCharge());
    ASSERT_TRUE(checkOrder(&scores, 1.0, 2.0, 4.0, 3.0, 5.0));

    // scan values are assigned [ 2, 2, 1, 1, 0 ]
    // labels are assigned [ -1, +1, -1, +1, -1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription();
        pPSM->scan = 2 - i / 2;
        scores.push_back(ScoreHolder(1.0 + i, (i % 2 ? LabelType::TARGET : LabelType::DECOY), pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanLabel());
    ASSERT_TRUE(checkOrder(&scores, 5.0, 4.0, 3.0, 2.0, 1.0));

    // scan values are assigned [ 2, 2, 1, 1, 0 ]
    // labels are assigned [ -1, +1, -1, +1, -1 ]
    // specFileNr values are assigned [ 0, 0, 1, 1, 2]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription();
        pPSM->specFileNr = i / 2;
        pPSM->scan = 2 - i / 2;
        scores.push_back(ScoreHolder(1.0 + i, (i % 2 ? LabelType::TARGET : LabelType::DECOY), pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanLabel());
    ASSERT_TRUE(checkOrder(&scores, 2.0, 1.0, 4.0, 3.0, 5.0));
}

// Verify the ScoreHolder's uniqueness functions.
TEST_F(ScoreHolderTest, CheckUniquenessFilter)
{
    std::vector<ScoreHolder> scores;

    // use unique peptides but repeating scan values
    // peptides are assigned in alphabetic order
    // scan values are assigned [ 0, 0, 1, 1, 2 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = newPsm(psmNames[i]);
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    ASSERT_EQ(5, scores.size());
    scores.erase(std::unique(scores.begin(), scores.end(), UniqueScanLabel()),
                 scores.end());
    ASSERT_EQ(3, scores.size());

    // use unique peptides but repeating specFileNr and scan values
    // peptides are assigned in alphabetic order
    // specFileNr values are assigned [ 0, 0, 0, 1, 1]
    // scan values are assigned [ 0, 0, 1, 1, 2 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = newPsm(psmNames[i]);
        pPSM->specFileNr = i / 3;
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    ASSERT_EQ(5, scores.size());
    scores.erase(std::unique(scores.begin(), scores.end(), UniqueScanLabel()),
                 scores.end());
    ASSERT_EQ(4, scores.size());
}


class ScoresTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
  private:
    int origVerbose;
};

void ScoresTest::SetUp()
{
    origVerbose = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(0);
}

void ScoresTest::TearDown()
{
    Globals::getInstance()->setVerbose(origVerbose);
}

// Test the basic populateWithPsms() method.
TEST_F(ScoresTest, CheckPopulating)
{
    Scores scores(true);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet *set1 = new DataSet();
    DataSet *set2 = new DataSet();
    set1->setLabel(LabelType::TARGET);
    set2->setLabel(LabelType::DECOY);
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *psm;
        psm = psmStore.add(PSMDescription());
        psmStore.storeStrings(*psm, "", psmNames[i], {});
        psm->scan = 100 + i;
        set1->registerPsm(psm);
        psm = psmStore.add(PSMDescription());
        psmStore.storeStrings(*psm, "", psmNames[i], {});
        psm->scan = 105 + i;
        set2->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set1));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set2));

    scores.populateWithPSMs(setHandler);
    ASSERT_EQ(10, scores.size());

    unsigned scanValue = 100;
    for (std::vector<ScoreHolder>::const_iterator it = scores.begin() ;
            it != scores.end() ;
            it++, scanValue++) {
        ASSERT_EQ(scanValue, it->pPSM->scan);
    }
}

// Test that populateWithPsms() refuses empty sets.
TEST_F(ScoresTest, CheckPopulatingEmpty)
{
    Scores scores(true);
    SetHandler setHandler(0);
    EXPECT_THROW(scores.populateWithPSMs(setHandler), MyException);
    DataSet *set1 = new DataSet();
    DataSet *set2 = new DataSet();
    set1->setLabel(LabelType::TARGET);
    set2->setLabel(LabelType::DECOY);
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set1));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set2));
    EXPECT_THROW(scores.populateWithPSMs(setHandler), MyException);
}

// createXvalSetsBySpectrum with N folds puts every PSM in the test set of
// exactly one fold and in the training sets of the other N - 1 folds, keeps
// the PSMs of a spectrum in one fold and gives folds of about equal size.
// With one fold, the training set and the test set both hold all PSMs.
TEST_F(ScoresTest, CreateXvalSetsBySpectrumWithNFolds)
{
    // 2 files with 37 spectra each, with 1 to 4 PSMs per spectrum and the
    // same scan numbers in both files
    std::vector<PSMDescription> psms;
    psms.reserve(400);
    for (unsigned int file = 0; file < 2; ++file) {
        for (unsigned int scan = 0; scan < 37; ++scan) {
            for (unsigned int k = 0; k <= (scan * 7 + file) % 4; ++k) {
                psms.emplace_back();
                psms.back().specFileNr = file;
                psms.back().scan = 1000 + scan;
            }
        }
    }
    const std::size_t maxPsmsPerSpectrum = 4;
    ASSERT_LT(psms.size(), 400u);  // no reallocation
    auto spectrumOf = [](const PSMDescription* psm) {
        return std::make_pair(psm->specFileNr, psm->scan);
    };

    for (unsigned int numFolds : {1u, 2u, 3u, 5u, 10u}) {
        Scores scores(false);
        for (std::size_t i = 0; i < psms.size(); ++i) {
            scores.addScoreHolder(ScoreHolder(
                0.0, i % 3 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
        }
        scores.recalculateSizes();
        ASSERT_EQ(74u, scores.getNumSpectra());

        std::vector<Scores> train, test;
        scores.createXvalSetsBySpectrum(train, test, numFolds);
        ASSERT_EQ(numFolds, train.size());
        ASSERT_EQ(numFolds, test.size());

        std::map<const PSMDescription*, unsigned int> testFold, numTrain;
        std::map<std::pair<unsigned int, unsigned int>, unsigned int>
            spectrumFold;
        std::size_t numTest = 0;
        for (unsigned int fold = 0; fold < numFolds; ++fold) {
            EXPECT_EQ(test[fold].size(),
                      static_cast<unsigned int>(
                          std::distance(test[fold].begin(), test[fold].end())));
            EXPECT_LE(test[fold].size(),
                      (psms.size() + numFolds - 1) / numFolds +
                          maxPsmsPerSpectrum - 1)
                << numFolds << " folds";
            EXPECT_GT(test[fold].size(), 0u) << numFolds << " folds";
            for (const ScoreHolder& sh : test[fold]) {
                ++numTest;
                EXPECT_EQ(0u, testFold.count(sh.pPSM));  // in one test set
                testFold[sh.pPSM] = fold;
                auto spectrum = spectrumOf(sh.pPSM);
                if (spectrumFold.count(spectrum)) {
                    EXPECT_EQ(spectrumFold[spectrum], fold);
                } else {
                    spectrumFold[spectrum] = fold;
                }
            }
        }
        EXPECT_EQ(psms.size(), numTest);
        EXPECT_EQ(74u, spectrumFold.size());
        for (unsigned int fold = 0; fold < numFolds; ++fold) {
            for (const ScoreHolder& sh : train[fold]) {
                ++numTrain[sh.pPSM];
                // a PSM is not trained on in the fold it is tested in, unless
                // there is only one fold
                if (numFolds > 1) {
                    EXPECT_NE(testFold[sh.pPSM], fold);
                }
            }
        }
        EXPECT_EQ(psms.size(), numTrain.size());
        for (const auto& psmCount : numTrain) {
            EXPECT_EQ(std::max(1u, numFolds - 1), psmCount.second);
        }
        if (numFolds == 1) {
            EXPECT_EQ(psms.size(), train[0].size());
            EXPECT_EQ(psms.size(), test[0].size());
            EXPECT_EQ(scores.posSize(), train[0].posSize());
            EXPECT_EQ(scores.negSize(), test[0].negSize());
        }
    }
}

// normalizeScores rescales the scores of each cross-validation bin to a
// common scale before the bins are merged. When the median decoy does not
// score below the score at the FDR threshold, which happens for tiny inputs,
// this is impossible and an error, unless the bin is the only one: its raw
// SVM scores and weights are then kept. A bin that can be rescaled is
// rescaled in the same way with and without singleBin.
TEST_F(ScoresTest, NormalizeScoresKeepsRawScoresOfASingleBin)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    const double fdr = 0.01;
    const std::vector<double> rawWeights = {0.5, -0.25, 1.0};
    std::vector<PSMDescription> psms(20);
    for (std::size_t i = 0; i < psms.size(); ++i) {
        psms[i].scan = static_cast<unsigned int>(i);
    }
    // a decoy above three targets, sorted by score; with a single decoy
    // there is no median decoy below the top score
    const double unscalableScores[] = {4.0, 3.0, 2.0, 1.0};
    auto unscalable = [&psms, &unscalableScores]() {
        Scores scores(false);
        for (std::size_t i = 0; i < 4; ++i) {
            scores.addScoreHolder(ScoreHolder(unscalableScores[i],
                i == 0 ? LabelType::DECOY : LabelType::TARGET, &psms[i]));
        }
        scores.recalculateSizes();
        return scores;
    };
    auto expectUnchanged = [&unscalableScores](const Scores& scores) {
        ASSERT_EQ(4u, scores.size());
        std::size_t i = 0;
        for (const ScoreHolder& sh : scores) {
            EXPECT_EQ(unscalableScores[i++], sh.score);
        }
    };

    {
        SCOPED_TRACE("normalizeScores");
        Scores scores = unscalable();
        scores.calcQvals(fdr);
        std::vector<double> w = rawWeights;
        EXPECT_THROW(scores.normalizeScores(fdr, w), MyException);
        EXPECT_NO_THROW(scores.normalizeScores(fdr, w, true));
        EXPECT_EQ(rawWeights, w);
        expectUnchanged(scores);
    }
    {
        SCOPED_TRACE("merge of one bin");
        std::vector<Scores> bins(1, unscalable());
        std::vector<std::vector<double> > w(1, rawWeights);
        Scores merged(false);
        EXPECT_NO_THROW(merged.merge(bins, fdr, false, w));
        EXPECT_EQ(rawWeights, w[0]);
        expectUnchanged(merged);
    }
    {
        SCOPED_TRACE("merge of two bins");
        std::vector<Scores> bins(2, unscalable());
        std::vector<std::vector<double> > w(2, rawWeights);
        Scores merged(false);
        EXPECT_THROW(merged.merge(bins, fdr, false, w), MyException);
    }
    {
        SCOPED_TRACE("rescalable bin");
        // 10 targets above 10 decoys
        Scores scores[2] = {Scores(false), Scores(false)};
        std::vector<double> w[2] = {rawWeights, rawWeights};
        for (int k = 0; k < 2; ++k) {
            for (std::size_t i = 0; i < 20; ++i) {
                scores[k].addScoreHolder(ScoreHolder(
                    20.0 - static_cast<double>(i),
                    i < 10 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
            }
            scores[k].recalculateSizes();
            scores[k].calcQvals(fdr);
            EXPECT_NO_THROW(scores[k].normalizeScores(fdr, w[k], k == 1));
        }
        EXPECT_NE(rawWeights, w[0]);
        EXPECT_EQ(w[0], w[1]);
        ASSERT_EQ(20u, scores[0].size());
        auto it0 = scores[0].begin();
        auto it1 = scores[1].begin();
        EXPECT_NE(20.0, it0->score);
        for (; it0 != scores[0].end(); ++it0, ++it1) {
            EXPECT_EQ(it0->score, it1->score);
        }
    }
    FeatureNames::resetNumFeatures();
}

// countTargetsBelowFdr must report the same number of targets below the FDR
// threshold as calcScoresAndQvals, also with many tied scores, with and
// without the mix-max correction.
TEST(ScoresCountTest, CountTargetsBelowFdrMatchesCalcScoresAndQvals)
{
    const size_t numFeatures = 2;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(numFeatures);
    const int numPsms = 2000;
    std::vector<double> features(numPsms * numFeatures);
    std::vector<PSMDescription*> psms;
    for (bool usePi0 : {false, true}) {
        Scores scores(usePi0);
        unsigned int seed = 12345u;
        for (int i = 0; i < numPsms; ++i) {
            seed = seed * 1103515245u + 12345u;
            bool isTarget = (seed >> 16) % 2 == 0;
            double* row = &features[i * numFeatures];
            // coarse values, so that many PSMs get exactly the same score
            row[0] = static_cast<double>((seed >> 8) % 7) + (isTarget ? 3 : 0);
            row[1] = static_cast<double>((seed >> 4) % 3);
            PSMDescription* psm = new PSMDescription();
            psm->features = row;
            psm->scan = i / 2;
            psms.push_back(psm);
            scores.addScoreHolder(ScoreHolder(0.0,
                isTarget ? LabelType::TARGET : LabelType::DECOY, psm));
        }
        scores.recalculateSizes();
        std::vector<double> w = {1.0, -0.5, 0.25};
        scores.calcScores(w);
        scores.postMergeStep();  // sets pi0 when usePi0 is true
        if (usePi0) {
            ASSERT_LT(scores.getPi0(), 1.0);
        }
        Scores::FeatureRows rows;
        scores.getFeatureRows(rows);
        Scores::CountBuffers buffers;
        // countTargetsBelowFdr needs the q values of the decoys as well;
        // calcScoresAndQvals sets this only after the first call below
        PosteriorEstimator::setNegative(true);
        for (double fdr : {0.01, 0.05, 0.2}) {
            for (bool skipDecoysPlusOne : {false, true}) {
                int count = scores.countTargetsBelowFdr(
                    rows, w, fdr, skipDecoysPlusOne, buffers);
                EXPECT_EQ(scores.calcScoresAndQvals(w, fdr, skipDecoysPlusOne),
                          count);
                if (fdr > 0.1) {
                    EXPECT_GT(count, 0);  // the test needs accepted targets
                }
            }
        }
    }
    for (PSMDescription* psm : psms) {
        delete psm;
    }
    FeatureNames::resetNumFeatures();
}

// Both PEP pathways of calcPep, from target q values (default) and from the
// target-decoy labels (tdc), with I-spline and PAVA regression, give PEPs in
// [0, 1] that do not decrease down the list, and whose mean over the targets
// above a q value threshold is close to that threshold.
TEST(ScoresPepTest, BothPathwaysGiveMonotoneCalibratedPeps)
{
    const int numPsms = 6000;
    std::vector<PSMDescription> psms(numPsms);
    unsigned int seed = 777u;
    auto uniform = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) + 0.5) / 16777216.0;
    };
    auto normal = [&uniform]() {
        return std::sqrt(-2.0 * std::log(uniform())) *
               std::cos(6.283185307179586 * uniform());
    };
    for (bool tdc : {false, true}) {
        for (bool pava : {false, true}) {
            Scores scores(false);
            seed = 777u;
            for (int i = 0; i < numPsms; ++i) {
                psms[i].scan = static_cast<unsigned int>(i);
                const bool isTarget = uniform() < 0.5;
                const bool correct = isTarget && uniform() < 0.4;
                scores.addScoreHolder(ScoreHolder(normal() + (correct ? 3.0 : 0.0),
                    isTarget ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
            }
            scores.recalculateSizes();
            // calcQvals and calcPep expect the PSMs sorted by score
            std::sort(scores.begin(), scores.end(), std::greater<ScoreHolder>());
            scores.calcQvals(0.01);
            scores.calcPep(false, true, pava, tdc);
            double lastTargetPep = 0.0, sum = 0.0;
            int n = 0;
            for (const ScoreHolder& sh : scores) {
                ASSERT_GE(sh.pep, 0.0);
                ASSERT_LE(sh.pep, 1.0);
                if (!sh.isTarget()) continue;
                EXPECT_GE(sh.pep, lastTargetPep - 1e-12);
                lastTargetPep = sh.pep;
                if (sh.q <= 0.05) {
                    sum += sh.pep;
                    ++n;
                }
            }
            ASSERT_GT(n, 100);
            EXPECT_GT(sum / n, 0.05 / 2.0) << "tdc=" << tdc << " pava=" << pava;
            EXPECT_LT(sum / n, 0.05 * 2.0) << "tdc=" << tdc << " pava=" << pava;
        }
    }
}

namespace {

// The number of targets with a q value below fdr as calcQvals counts it
// with pi0 = 1: sort all PSMs, best first, and count with getQValues.
int referenceTargetCount(std::vector<std::pair<double, bool> > psms,
                         bool lowBest, double fdr, bool skipDecoysPlusOne,
                         double nullTargetWinProb) {
    std::sort(psms.begin(), psms.end(),
              [lowBest](const std::pair<double, bool>& a,
                        const std::pair<double, bool>& b) {
                  return lowBest ? a.first < b.first : a.first > b.first;
              });
    std::vector<double> q;
    PosteriorEstimator::setNegative(true);
    PosteriorEstimator::getQValues(1.0, psms, q, skipDecoysPlusOne,
                                   nullTargetWinProb);
    int count = 0;
    for (std::size_t i = 0; i < q.size(); ++i) {
        if (q[i] < fdr && psms[i].second) ++count;
    }
    return count;
}

// Rows with the labels of psms; countTargetsBySelection reads only the
// labels from them, so they point to placeholder storage.
Scores::FeatureRows rowsFor(const std::vector<std::pair<double, bool> >& psms,
                            std::vector<double>& storage) {
    storage.assign(psms.size() + 1, 0.0);
    Scores::FeatureRows rows(psms.size());
    for (std::size_t i = 0; i < psms.size(); ++i) {
        rows[i] = Scores::makeFeatureRow(&storage[i], psms[i].second);
    }
    return rows;
}

// Compares both counts of countTargetsBySelection, and each count alone,
// with the full sort.
void expectSelectionCountsMatch(
        const std::vector<std::pair<double, bool> >& psms, double fdr,
        bool skipDecoysPlusOne, double nullTargetWinProb,
        Scores::CountBuffers& buffers, const std::string& what) {
    std::vector<double> storage;
    Scores::FeatureRows rows = rowsFor(psms, storage);
    buffers.values.resize(psms.size());
    for (std::size_t i = 0; i < psms.size(); ++i) {
        buffers.values[i] = psms[i].first;
    }
    const std::string where = what + " fdr=" + std::to_string(fdr) +
        " skip=" + std::to_string(skipDecoysPlusOne) + " p=" +
        std::to_string(nullTargetWinProb);
    int high = -2, low = -2;
    ASSERT_TRUE(Scores::countTargetsBySelection(
        rows, buffers, fdr, skipDecoysPlusOne, nullTargetWinProb, &high,
        &low)) << where;
    ASSERT_EQ(referenceTargetCount(psms, false, fdr, skipDecoysPlusOne,
                                   nullTargetWinProb), high)
        << where << " high first";
    ASSERT_EQ(referenceTargetCount(psms, true, fdr, skipDecoysPlusOne,
                                   nullTargetWinProb), low)
        << where << " low first";
    int single = -2;
    ASSERT_TRUE(Scores::countTargetsBySelection(
        rows, buffers, fdr, skipDecoysPlusOne, nullTargetWinProb, &single,
        nullptr)) << where;
    ASSERT_EQ(high, single) << where << " high first only";
    ASSERT_TRUE(Scores::countTargetsBySelection(
        rows, buffers, fdr, skipDecoysPlusOne, nullTargetWinProb, nullptr,
        &single)) << where;
    ASSERT_EQ(low, single) << where << " low first only";
}

}  // namespace

// Every label sequence of up to 10 PSMs, with distinct values, all values
// tied, blocks of 2 and 3 tied values and two alternating values, for FDR
// thresholds that the estimates reach exactly, both decoy offsets and the
// null target win probabilities of normal training and of RESET.
TEST(ScoresCountTest, SelectionCountMatchesFullSortForAllSmallInputs)
{
    OneThread oneThread;
    Scores::CountBuffers buffers;
    const double fdrs[] = {0.0, 0.1, 0.2, 0.25, 1.0 / 3.0, 0.5, 2.0 / 3.0,
                           1.0, 1.5};
    const double nullTargetWinProbs[] = {0.5, 0.75, 2.0 / 3.0, 0.0};
    std::vector<std::pair<double, bool> > psms;
    for (int n = 0; n <= 10; ++n) {
        for (unsigned int mask = 0; mask < (1u << n); ++mask) {
            for (int pattern = 0; pattern < 5; ++pattern) {
                psms.resize(n);
                for (int i = 0; i < n; ++i) {
                    double value = 0.0;
                    switch (pattern) {
                        case 0: value = i; break;
                        case 1: value = 0.0; break;
                        case 2: value = i / 2; break;
                        case 3: value = -(i / 3); break;
                        case 4: value = i % 2; break;
                    }
                    psms[i] = std::make_pair(value, ((mask >> i) & 1u) != 0);
                }
                const std::string what = "n=" + std::to_string(n) +
                    " mask=" + std::to_string(mask) + " pattern=" +
                    std::to_string(pattern);
                for (double fdr : fdrs) {
                    for (bool skip : {false, true}) {
                        for (double p : nullTargetWinProbs) {
                            expectSelectionCountsMatch(psms, fdr, skip, p,
                                                       buffers, what);
                            if (HasFatalFailure()) return;
                        }
                    }
                }
            }
        }
    }
}

// Larger random inputs in random order: coarse values with many ties, a
// few signed zeros and infinities, and all-target and all-decoy sets.
TEST(ScoresCountTest, SelectionCountMatchesFullSortForRandomInputs)
{
    OneThread oneThread;
    Scores::CountBuffers buffers;
    unsigned int seed = 4242u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return seed >> 8;
    };
    std::vector<std::pair<double, bool> > psms;
    for (int round = 0; round < 60; ++round) {
        const int n = 1 + static_cast<int>(next() % 4000);
        const unsigned int targetPercent =
            round % 10 == 0 ? 100 : (round % 10 == 1 ? 0 : next() % 100);
        const unsigned int levels = 2 + next() % 200;
        psms.resize(n);
        for (int i = 0; i < n; ++i) {
            const bool isTarget = next() % 100 < targetPercent;
            double value = static_cast<double>(next() % levels) +
                           (isTarget ? static_cast<double>(next() % 20) : 0.0);
            if (next() % 500 == 0) value = (next() % 2) ? 0.0 : -0.0;
            if (next() % 1000 == 0) {
                value = (next() % 2) ? std::numeric_limits<double>::infinity()
                                     : -std::numeric_limits<double>::infinity();
            }
            psms[i] = std::make_pair(value, isTarget);
        }
        for (double fdr : {0.001, 0.01, 0.05, 0.3}) {
            for (bool skip : {false, true}) {
                for (double p : {0.5, 0.75}) {
                    expectSelectionCountsMatch(psms, fdr, skip, p, buffers,
                                               "round " + std::to_string(round));
                    if (HasFatalFailure()) return;
                }
            }
        }
    }
}

// NaN values, and null target win probabilities without a finite decoy
// factor, are left to the full sort.
TEST(ScoresCountTest, SelectionCountRejectsNanAndInvalidProbabilities)
{
    Scores::CountBuffers buffers;
    std::vector<std::pair<double, bool> > psms = {
        {1.0, true}, {0.5, false}, {2.0, true}, {0.1, false}};
    std::vector<double> storage;
    Scores::FeatureRows rows = rowsFor(psms, storage);
    buffers.values = {1.0, 0.5, 2.0, 0.1};
    int high = -2, low = -2;
    EXPECT_TRUE(Scores::countTargetsBySelection(rows, buffers, 0.5, true, 0.5,
                                                &high, &low));
    EXPECT_EQ(2, high);
    EXPECT_EQ(0, low);
    EXPECT_FALSE(Scores::countTargetsBySelection(rows, buffers, 0.5, true, 1.0,
                                                 &high, &low));
    EXPECT_FALSE(Scores::countTargetsBySelection(rows, buffers, 0.5, true, 1.5,
                                                 &high, &low));
    buffers.values[1] = std::nan("");
    EXPECT_FALSE(Scores::countTargetsBySelection(rows, buffers, 0.5, true, 0.5,
                                                 &high, &low));
}

// getInitDirection picks the same feature and direction, with the same
// count, as the full sort of every feature column in the order of the
// PSMs (the earlier implementation): with ties, an informative feature in
// each direction, constant and NaN columns (the NaN column takes the full
// sort) and the null target win probabilities of normal training and of
// RESET.
TEST(ScoresCountTest, InitDirectionMatchesFullSortOfEachFeature)
{
    OneThread oneThread;
    const int origVerbose = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(0);
    const size_t numFeatures = 5;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(numFeatures);
    unsigned int seed = 99u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return seed >> 8;
    };
    for (int round = 0; round < 12; ++round) {
        const int numPsms = 50 + static_cast<int>(next() % 3000);
        std::vector<double> features(numPsms * numFeatures);
        std::vector<PSMDescription> psms(numPsms);
        Scores scores(false);
        const bool withNan = round % 3 == 2;
        for (int i = 0; i < numPsms; ++i) {
            const bool isTarget = next() % 2 == 0;
            const bool correct = isTarget && next() % 3 == 0;
            double* row = &features[i * numFeatures];
            row[0] = static_cast<double>(next() % 50) + (correct ? 20 : 0);
            row[1] = -static_cast<double>(next() % 7) - (correct ? 3 : 0);
            row[2] = static_cast<double>(next() % 1000) / 7.0;
            row[3] = 1.0;
            row[4] = withNan && next() % 100 == 0
                         ? std::nan("")
                         : static_cast<double>(next() % 30) + (correct ? 5 : 0);
            psms[i].features = row;
            psms[i].scan = static_cast<unsigned int>(i);
            scores.addScoreHolder(ScoreHolder(0.0,
                isTarget ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
        }
        scores.recalculateSizes();
        const double nullTargetWinProb = round % 2 == 0 ? 0.5 : 0.75;
        scores.setNullTargetWinProb(nullTargetWinProb);
        for (double fdr : {0.01, 0.1, 0.3}) {
            // the earlier implementation, which sorts every column
            int bestPositives = -1, bestFeature = -1;
            bool lowBest = false;
            for (size_t f = 0; f < numFeatures; ++f) {
                std::vector<std::pair<double, bool> > combined;
                for (const ScoreHolder& sh : scores) {
                    combined.push_back(std::make_pair(sh.pPSM->features[f],
                                                      sh.isTarget()));
                }
                std::sort(combined.begin(), combined.end(),
                          [](const std::pair<double, bool>& a,
                             const std::pair<double, bool>& b) {
                              return a.first < b.first;
                          });
                for (int i = 0; i < 2; ++i) {
                    if (i == 1) std::reverse(combined.begin(), combined.end());
                    std::vector<double> q;
                    PosteriorEstimator::setNegative(true);
                    PosteriorEstimator::getQValues(1.0, combined, q, true,
                                                   nullTargetWinProb);
                    int positives = 0;
                    for (size_t k = 0; k < q.size(); ++k) {
                        if (q[k] < fdr && combined[k].second) ++positives;
                    }
                    if (positives > bestPositives) {
                        bestPositives = positives;
                        bestFeature = static_cast<int>(f);
                        lowBest = (i == 0);
                    }
                }
            }
            ASSERT_GT(bestPositives, 0);
            std::vector<double> direction(numFeatures + 1, 7.0);
            EXPECT_EQ(bestPositives, scores.getInitDirection(fdr, direction))
                << "round " << round << " fdr " << fdr;
            for (size_t f = 0; f < numFeatures; ++f) {
                const double expected = static_cast<int>(f) == bestFeature
                                            ? (lowBest ? -1.0 : 1.0)
                                            : 0.0;
                EXPECT_EQ(expected, direction[f])
                    << "round " << round << " fdr " << fdr << " feature " << f;
            }
        }
    }
    FeatureNames::resetNumFeatures();
    Globals::getInstance()->setVerbose(origVerbose);
}

// isSortedByScore is true exactly when a sort by decreasing score, ties
// broken by the total order of ScoreHolder, would keep the order; calcScores
// and merge skip that sort then.
TEST(ScoresSortTest, IsSortedByScoreMatchesTheOrderOfTheSort)
{
    std::vector<PSMDescription> psms(6);
    for (size_t i = 0; i < psms.size(); ++i) {
        psms[i].scan = static_cast<unsigned int>(10 - i);
    }
    auto makeScores = [&psms](const std::vector<double>& values) {
        Scores scores(false);
        for (size_t i = 0; i < values.size(); ++i) {
            scores.addScoreHolder(ScoreHolder(values[i],
                i % 2 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
        }
        return scores;
    };
    EXPECT_TRUE(makeScores({}).isSortedByScore());
    EXPECT_TRUE(makeScores({1.0}).isSortedByScore());
    EXPECT_TRUE(makeScores({3.0, 2.0, 2.0, 1.0}).isSortedByScore());
    EXPECT_FALSE(makeScores({3.0, 1.0, 2.0}).isSortedByScore());
    EXPECT_FALSE(makeScores({1.0, 2.0}).isSortedByScore());
    EXPECT_FALSE(makeScores({3.0, std::nan(""), 1.0}).isSortedByScore());
    EXPECT_FALSE(makeScores({3.0, 2.0, std::nan("")}).isSortedByScore());
    // ties are ordered by decreasing scan number: psms[i] has scan 10 - i
    Scores tied = makeScores({2.0, 2.0, 2.0});
    EXPECT_TRUE(tied.isSortedByScore());
    std::swap(*tied.begin(), *(tied.begin() + 2));
    EXPECT_FALSE(tied.isSortedByScore());
    std::sort(tied.begin(), tied.end(), std::greater<ScoreHolder>());
    EXPECT_TRUE(tied.isSortedByScore());
    // a sorted set stays as it is when scored again with the same weights
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(1);
    std::vector<double> features = {0.5, 0.25, 0.25, 1.0, -1.0, 0.0};
    Scores scored(false);
    for (size_t i = 0; i < psms.size(); ++i) {
        psms[i].features = &features[i];
        scored.addScoreHolder(ScoreHolder(0.0, LabelType::TARGET, &psms[i]));
    }
    std::vector<double> w = {2.0, 0.5};
    scored.calcScores(w);
    ASSERT_TRUE(scored.isSortedByScore());
    std::vector<const PSMDescription*> order;
    for (const ScoreHolder& sh : scored) order.push_back(sh.pPSM);
    scored.calcScores(w);
    size_t k = 0;
    for (const ScoreHolder& sh : scored) EXPECT_EQ(order[k++], sh.pPSM);
    FeatureNames::resetNumFeatures();
}

namespace {

// Scores with access to the protected members that the q value and p value
// tests set or compare.
class ExposedScores : public Scores {
 public:
  explicit ExposedScores(bool usePi0) : Scores(usePi0) {}
  void setPi0(double pi0) { pi0_ = pi0; }
  using Scores::getScoreLabelPairs;
  using Scores::getTargetPValues;
};

// calcQvals must give bit for bit the q values and the count of
// PosteriorEstimator::getQValues on the (score, isTarget) pairs, including
// the out_of_range of its mix-max lookup with skipDecoysPlusOne when the
// best-scoring group with decoys has a single decoy. Returns whether
// getQValues threw.
bool expectQvalsAsGetQValues(ExposedScores& scores, double pi0,
                             double nullTargetWinProb, bool skipDecoysPlusOne,
                             double fdr, const std::string& what) {
  scores.setPi0(pi0);
  scores.setNullTargetWinProb(nullTargetWinProb);
  std::vector<std::pair<double, bool> > combined;
  scores.getScoreLabelPairs(combined);
  std::vector<double> expected;
  PosteriorEstimator::setNegative(true);
  bool expectThrow = false;
  try {
    PosteriorEstimator::getQValues(pi0, combined, expected, skipDecoysPlusOne,
                                   nullTargetWinProb);
  } catch (const std::out_of_range&) {
    expectThrow = true;
  }
  if (expectThrow) {
    EXPECT_THROW(scores.calcQvals(fdr, skipDecoysPlusOne), std::out_of_range)
        << what;
    return true;
  }
  const int numPos = scores.calcQvals(fdr, skipDecoysPlusOne);
  EXPECT_EQ(expected.size(), scores.size()) << what;
  int expectedPos = 0;
  size_t ix = 0;
  for (const ScoreHolder& sh : scores) {
    if (ix >= expected.size() || expected[ix] != sh.q) {
      ADD_FAILURE() << what << ", position " << ix;
      return false;
    }
    if (expected[ix] < fdr && combined[ix].second) ++expectedPos;
    ++ix;
  }
  EXPECT_EQ(expectedPos, numPos) << what;
  return false;
}

void expectPValuesAsGetPValues(ExposedScores& scores, const std::string& what) {
  std::vector<std::pair<double, bool> > combined;
  scores.getScoreLabelPairs(combined);
  std::vector<double> expected, pvals;
  PosteriorEstimator::getPValues(combined, expected);
  scores.getTargetPValues(pvals);
  ASSERT_EQ(expected.size(), pvals.size()) << what;
  for (size_t ix = 0; ix < pvals.size(); ++ix) {
    ASSERT_EQ(expected[ix], pvals[ix]) << what << ", position " << ix;
  }
}

const double kPi0s[] = {1.0, 0.7, 0.3};
const double kNullTargetWinProbs[] = {0.5, 2.0 / 3.0, 0.75};

}  // namespace

// Every target/decoy pattern of up to 8 PSMs, with distinct scores, all
// scores tied and blocks of tied scores.
TEST(ScoresQvalueTest, AllSmallPatternsEqualGetQValues)
{
    std::vector<PSMDescription> psms(8);
    int numThrows = 0;
    for (int n = 0; n <= 8; ++n) {
        for (int pattern = 0; pattern < (1 << n); ++pattern) {
            for (int block : {1, 2, 3, 8}) {
                ExposedScores scores(true);
                for (int i = 0; i < n; ++i) {
                    const bool isTarget = (pattern >> i) & 1;
                    scores.addScoreHolder(ScoreHolder(
                        -static_cast<double>(i / block),
                        isTarget ? LabelType::TARGET : LabelType::DECOY,
                        &psms[i]));
                }
                scores.recalculateSizes();
                std::string what = "n=" + std::to_string(n) + " pattern=" +
                                   std::to_string(pattern) +
                                   " block=" + std::to_string(block);
                expectPValuesAsGetPValues(scores, what);
                for (double pi0 : kPi0s) {
                    for (double nullTargetWinProb : kNullTargetWinProbs) {
                        for (bool skip : {false, true}) {
                            numThrows += expectQvalsAsGetQValues(
                                scores, pi0, nullTargetWinProb, skip, 0.3,
                                what + " pi0=" + std::to_string(pi0) +
                                    " ntwp=" + std::to_string(nullTargetWinProb) +
                                    " skip=" + std::to_string(skip));
                        }
                    }
                }
            }
        }
    }
    EXPECT_GT(numThrows, 0);  // the out_of_range case was among them
}

// Larger random sets with many ties and pseudo targets (RESET), including
// pi0 and nullTargetWinProb as the training sets and postMergeStep set them.
TEST(ScoresQvalueTest, RandomSetsEqualGetQValues)
{
    const int numPsms = 20000;
    std::vector<PSMDescription> psms(numPsms);
    unsigned int seed = 2024u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return seed >> 8;
    };
    for (int coarseness : {3, 40, 1000000}) {
        ExposedScores scores(true);
        for (int i = 0; i < numPsms; ++i) {
            const unsigned int r = next();
            const LabelType label = r % 5 == 0   ? LabelType::PSEUDO_TARGET
                                    : r % 5 < 3 ? LabelType::TARGET
                                                : LabelType::DECOY;
            const bool good = label == LabelType::TARGET && next() % 3 == 0;
            const double score =
                static_cast<double>(next() % coarseness) / coarseness +
                (good ? 0.5 : 0.0);
            psms[i].scan = static_cast<unsigned int>(i);
            scores.addScoreHolder(ScoreHolder(score, label, &psms[i]));
        }
        scores.recalculateSizes();
        std::sort(scores.begin(), scores.end(), std::greater<ScoreHolder>());
        std::string what = "coarseness=" + std::to_string(coarseness);
        expectPValuesAsGetPValues(scores, what);
        for (double pi0 : kPi0s) {
            for (double nullTargetWinProb : kNullTargetWinProbs) {
                for (bool skip : {false, true}) {
                    for (double fdr : {0.01, 0.1}) {
                        expectQvalsAsGetQValues(
                            scores, pi0, nullTargetWinProb, skip, fdr,
                            what + " pi0=" + std::to_string(pi0) +
                                " ntwp=" + std::to_string(nullTargetWinProb) +
                                " skip=" + std::to_string(skip) +
                                " fdr=" + std::to_string(fdr));
                    }
                }
            }
        }
        // the pi0 that postMergeStep estimates from these p values
        scores.postMergeStep();
        ASSERT_LT(scores.getPi0(), 1.0);
        expectQvalsAsGetQValues(scores, scores.getPi0(), 0.5, false, 0.01,
                                what + " estimated pi0");
    }
}
