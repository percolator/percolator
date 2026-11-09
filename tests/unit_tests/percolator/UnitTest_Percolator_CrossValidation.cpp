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

#include <gtest/gtest.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>
#include "SetHandler.h"
#include "CrossValidation.h"
#include "Globals.h"
#include "PseudoRandom.h"
#include "NoNormalizer.h"
#include <cstdio>
#include <fstream>

/* A subclass of CrossValidation that gives us access to some
 * protected fields.
 */
class CrossValidationEx : public CrossValidation {
public:
    CrossValidationEx(bool quickValidation,
                      bool reportPerformanceEachIteration, double testFdr,
                      double selectionFdr, double initialSelectionFdr,
                      double selectedCpos, double selectedCneg,
                      unsigned int niter, bool usePi0,
                      unsigned int nestedXvalBins, bool trainBestPositive,
                      unsigned int numThreads, bool skipNormalizeScores, 
                      double decoyFractionTraining, unsigned int numFolds)
        : CrossValidation(quickValidation, reportPerformanceEachIteration,
                          testFdr, selectionFdr, initialSelectionFdr,
                          selectedCpos, selectedCneg, niter, usePi0,
                          nestedXvalBins, trainBestPositive, numThreads,
                          skipNormalizeScores, decoyFractionTraining, numFolds) { }
    virtual ~CrossValidationEx() { }
    int doStepEx(Normalizer* pNorm, double selectionFdr) {
        return doStep(pNorm, selectionFdr);
    }
    std::vector< std::vector<double> > const& weights(void) const {
        return weights_;
    }
    std::vector<Scores> const& testBins(void) const {
        return testScores_;
    }
    std::vector<Scores> const& trainBins(void) const {
        return trainScores_;
    }
};

/*
 * Unit tests for the CrossValidation class.
 *
 * Note: These tests take over 250 ms to complete (on my computer).
 * They are okay here as long as these are the only slow tests in the
 * unit test suite. But at some point we will probably want to move
 * this class into a separate suite of integration tests.
 */

class CrossValidationTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
    int origVerbose;
};

void CrossValidationTest::SetUp()
{
    origVerbose = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(0);
}

void CrossValidationTest::TearDown()
{
    Globals::getInstance()->setVerbose(origVerbose);
}

// With numFolds folds each PSM is in numFolds - 1 training sets, and the
// estimated number of positives is their sum over the training sets divided
// by numFolds - 1; with one fold the only training set holds all PSMs. The
// estimate is therefore about the number of positives in the whole data set
// (N) for any number of folds.
TEST_F(CrossValidationTest, doStepTest)
{
  // The folds are drawn with the process-wide generator, which other tests
  // also draw from; start from the default seed so that the folds, and the
  // estimates checked below, do not depend on which tests ran before.
  PseudoRandom::setSeed(1u);
  // the default 3 folds first, then other numbers of folds; with 2 folds the
  // training sets hold only about N/2 targets above the first decoy, too few
  // to reach q < testFdr in some of them, and the estimate is about N/2
  for (unsigned int numFolds : {3u, 1u, 5u, 10u}) {
    SCOPED_TRACE(numFolds);
    // Note that we use an elevated testFdr (0.02 instead of 0.01),
    // to compensate for the small size of the data sets.
    int const N = 100;
    double const testFdr = 0.02;

    CrossValidationEx *crossValidation =
            new CrossValidationEx(false,   // quickValidation
                                  false,   // reportEachIteration
                                  testFdr, // testFdr
                                  0.01,    // selectionFdr
                                  0.01,    // initialSectionFdr
                                  0.0,     // selectedCpos
                                  0.0,     // selectedCneg
                                  10,      // nIter
                                  true,    // usePi0
                                  1,       // nestedXvalBins
                                  false,   // trainBestPositive
                                  1,       // numThreads
                                  false,   // skipNormalizeScores
                                  1.0,     // decoyFractionTraining
                                  numFolds);

    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    SetHandler setHandler(0);
    int scanNumber = 1;

    // Our data set has two features. One flips between 0 and 1
    // without regard for label. The other is consistently increasing,
    // but with the targets' values slightly less than the decoys'
    // values.

    PsmStore psmStore;
    // a PSMDescription does not own its features
    std::vector<double> featureRows(3 * N * 2);
    double *nextRow = featureRows.data();
    DataSet *targets = new DataSet();
    DataSet *decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    for (int i = 0 ; i < 2 * N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = nextRow;
        nextRow += 2;
        psm->features[0] = static_cast<double>(i) / N - 1E-6;
        psm->features[1] = static_cast<double>(i % 2);
        psm->scan = scanNumber++;
        targets->registerPsm(psm);
    }
    for (int i = 0 ; i < N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = nextRow;
        nextRow += 2;
        psm->features[0] = static_cast<double>(i) / N;
        psm->features[1] = static_cast<double>(i % 2);
        psm->scan = scanNumber++;
        decoys->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));

    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);
    Normalizer* pNorm = norm.get();

    Scores scores(true);
    scores.populateWithPSMs(setHandler);

    // The pre-iteration setup should select the first feature to
    // be weighted.

    SanityCheck check;
    int numInit = crossValidation->preIterationSetup(scores, check, pNorm);
    std::vector< std::vector<double> > const& w = crossValidation->weights();
    ASSERT_EQ(numFolds, w.size());
    for (unsigned int i = 0 ; i < numFolds ; ++i) {
        EXPECT_EQ(1.0, w[i][0]);
        EXPECT_EQ(0.0, w[i][1]);
    }

    // One step of the training algorithm should find at least N positives.
    // The number can be a bit higher because some targets with score < N are
    // counted as positives. This can be either because the highest scoring decoy 
    // in the cross validation fold has a lower score, or because it narrowly passes 
    // the testFDR threshold.
    // The estimate is truncated, and the trained direction can move single
    // targets below a decoy in some training sets, so it can be slightly
    // below N for other numbers of folds.
    int estimatedNumPositives = crossValidation->doStepEx(pNorm, 0.01);
    if (numFolds == 3u) {
      EXPECT_LE(N, estimatedNumPositives);
    } else {
      EXPECT_LE(N * (1.0 - testFdr), estimatedNumPositives);
    }
    EXPECT_GE(N * (1.0 + testFdr), estimatedNumPositives);

    delete crossValidation;
  }
  FeatureNames::resetNumFeatures();
}

// The direction check and the data of trainAndMerge
struct DirectionCheckOptions {
    // --direction-check
    SanityCheck::DirectionCheck check = SanityCheck::DirectionCheck::BINS;
    bool overrule = false;  // --override
    double testFdr = 0.05;  // a tiny one makes the check fail
    // the merge after the second chance kept the trained SVMs
    MergeNormalization secondChanceScheme = MergeNormalization::TRAIN_DECOYS;
    bool signal = true;     // without, no target is correct
    bool separate = false;  // pi0 and q values of separate searches (mix-max)
    // the sum of the two features separates every target from every decoy,
    // while neither feature does alone
    bool perfect = false;
    // --static: the scores of the bins are not rescaled
    bool skipNormalizeScores = false;
    // --init-weights: a weights file whose raw weights are the initial
    // direction of every fold
    std::string initWeightFN;
};

// What trainAndMerge observed besides the scores and weights
struct DirectionCheckResult {
    bool secondChanceKept = false;
    MergeNormalization usedScheme = MergeNormalization::TEST;
    // validateDirection, run on copies before postIterationProcessing,
    // left the random number generator as it found it
    bool checkKeptRandomState = false;
    // the state after postIterationProcessing, whose merge draws numbers
    // for the pi0 of separate searches
    uint64_t randomStateAfterMerge = 0u;
    // estimating the pi0 of the merged list once more moved the generator
    bool pi0EstimateDraws = false;
    // the (input line, fold) list of -N: the test PSMs that keep their
    // test fold's score
    std::size_t numTestFolds = 0u;
};

// Trains and merges with numFolds folds and the given merge scheme on 450
// targets and 450 decoys. Half of the targets score higher in both features,
// so that the SVM improves on the initial direction of a single feature.
// With checkOptions.perfect the targets get the features 2 + u and 2 - u,
// with u in [-4, 4), and the decoys features in [-1, 1).
// Returns the merged score of each scan and the final weights, and checks
// that the final weights of each fold give the merged scores of its test
// bin: each PSM was scored by the fold of its own test bin.
static void trainAndMerge(unsigned int numFolds, MergeNormalization scheme,
                          std::map<unsigned int, double>& scoreOf,
                          std::vector<std::vector<double> >& weights,
                          std::vector<double>& features,
                          const DirectionCheckOptions& checkOptions =
                              DirectionCheckOptions(),
                          DirectionCheckResult* result = nullptr) {
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    const int numPsms = 900;
    features.assign(2 * numPsms, 0.0);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet* targets = new DataSet();
    DataSet* decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    unsigned int seed = 2024u;
    for (int i = 0; i < numPsms; ++i) {
        seed = seed * 1103515245u + 12345u;
        const bool isTarget = i % 2 == 0;
        PSMDescription* psm = psmStore.add(PSMDescription());
        psm->features = &features[2 * i];
        const bool isCorrect = checkOptions.signal && isTarget && i % 4 == 0;
        psm->features[0] = static_cast<double>((seed >> 8) % 1000) / 250.0 +
                           (isCorrect ? 2.5 : 0.0);
        psm->features[1] = static_cast<double>((seed >> 4) % 1000) / 250.0 +
                           (isCorrect ? 2.5 : 0.0);
        if (checkOptions.perfect) {
            const double u = static_cast<double>((seed >> 8) % 1000) / 125.0;
            const double v = static_cast<double>((seed >> 4) % 1000) / 500.0;
            psm->features[0] = isTarget ? u - 2.0 : u / 4.0 - 1.0;
            psm->features[1] = isTarget ? 6.0 - u : v - 1.0;
        }
        psm->scan = static_cast<unsigned int>(1 + i);
        (isTarget ? targets : decoys)->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);

    PseudoRandom::setSeed(1u);
    CrossValidationEx crossValidation(false, false, checkOptions.testFdr,
                                      0.05, 0.05, 0.0, 0.0, 10,
                                      checkOptions.separate, 1, false, 1,
                                      checkOptions.skipNormalizeScores, 1.0,
                                      numFolds);
    crossValidation.setMergeNormalization(scheme);
    Scores scores(checkOptions.separate);
    scores.populateWithPSMs(setHandler);
    SanityCheck check;
    check.setDirectionCheck(checkOptions.check);
    check.setOverrule(checkOptions.overrule);
    check.setSecondChanceMergeScheme(checkOptions.secondChanceScheme);
    if (!checkOptions.initWeightFN.empty()) {
        check.setInitWeightFN(checkOptions.initWeightFN);
    }
    crossValidation.setRecordTestFolds(true);
    crossValidation.preIterationSetup(scores, check, norm.get());
    crossValidation.train(norm.get());
    // the bin of each PSM, before the merge empties the bins: the weights
    // give the merged scores, as each PSM was scored by the fold of its test
    // bin
    std::map<const PSMDescription*, std::size_t> binOf;
    for (std::size_t bin = 0; bin < numFolds; ++bin) {
        for (const ScoreHolder& sh : crossValidation.testBins()[bin]) {
            binOf[sh.pPSM] = bin;
        }
    }
    if (result) {
        SanityCheck checkCopy(check);
        std::vector<Scores> testCopy(crossValidation.testBins());
        std::vector<Scores> trainCopy(crossValidation.trainBins());
        std::vector<std::vector<double> > wCopy(crossValidation.weights());
        const uint64_t before = PseudoRandom::getSeed();
        checkCopy.validateDirection(testCopy, trainCopy, wCopy);
        result->checkKeptRandomState = PseudoRandom::getSeed() == before;
    }
    crossValidation.postIterationProcessing(scores, check);
    if (result) {
        result->secondChanceKept = check.secondChanceKept();
        result->usedScheme = crossValidation.usedMergeNormalization();
        result->randomStateAfterMerge = PseudoRandom::getSeed();
        std::vector<std::pair<unsigned int, unsigned int> > testFolds;
        crossValidation.getTestFolds(testFolds);
        result->numTestFolds = testFolds.size();
    }

    scoreOf.clear();
    for (const ScoreHolder& sh : scores) {
        scoreOf[sh.pPSM->scan] = sh.score;
    }
    weights = crossValidation.weights();
    ASSERT_EQ(static_cast<std::size_t>(numPsms), binOf.size());
    ASSERT_EQ(static_cast<unsigned int>(numPsms), scores.size());
    for (const ScoreHolder& sh : scores) {
        const double expected =
            scores.calcScore(sh.pPSM->features, weights[binOf[sh.pPSM]]);
        EXPECT_NEAR(expected, sh.score,
                    1e-12 * std::max(1.0, std::fabs(sh.score)))
            << "bin " << binOf[sh.pPSM];
    }
    if (result) {
        const uint64_t before = PseudoRandom::getSeed();
        scores.postMergeStep(false);
        result->pi0EstimateDraws = PseudoRandom::getSeed() != before;
    }
    FeatureNames::resetNumFeatures();
}

// trainAndMerge finds each test bin scored with its own fold's final
// weights. When the check passes, they are the trained ones. When it fails,
// here by a test FDR that no target reaches, they are, with --override, the
// trained ones, which were not trained on that bin, and otherwise the
// default vector of a single feature chosen on the fold's training set.
// Before, all bins were scored with the first fold's vector after a failed
// check; after a reset this failed even where the folds chose the same
// feature, since each fold's default vector kept its trained bias.
TEST_F(CrossValidationTest, FailedDirectionCheckScoresEachBinWithItsOwnVector)
{
    auto isSingleFeature = [](const std::vector<double>& w) {
        return (w[0] == 0.0) != (w[1] == 0.0);
    };
    std::map<unsigned int, double> scoreOf;
    std::vector<double> features;
    for (unsigned int numFolds : {3u, 5u}) {
        {
            SCOPED_TRACE(::testing::Message() << numFolds
                                              << " folds, check passes");
            std::vector<std::vector<double> > w;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                          features);
            ASSERT_EQ(numFolds, w.size());
            for (const std::vector<double>& wFold : w) {
                EXPECT_FALSE(isSingleFeature(wFold));  // trained, not reset
            }
        }
        for (bool overrule : {false, true}) {
            SCOPED_TRACE(::testing::Message() << numFolds
                                              << " folds, failed check, "
                                              << "override " << overrule);
            DirectionCheckOptions failing;
            failing.overrule = overrule;
            failing.testFdr = 1e-6;
            std::vector<std::vector<double> > w;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                          features, failing);
            ASSERT_EQ(numFolds, w.size());
            for (const std::vector<double>& wFold : w) {
                EXPECT_EQ(!overrule, isSingleFeature(wFold));
            }
            if (overrule) {
                // each fold keeps its own trained weights
                for (unsigned int fold = 1; fold < numFolds; ++fold) {
                    EXPECT_NE(w[0], w[fold]);
                }
            }
        }
    }
}

// With one fold the training set is the test bin, so the training-set
// anchors of --merge-normalization train are the test-bin anchors: the
// merged scores and the weights equal those of the default scheme exactly.
// With three folds every scheme gives scores that the final weights
// describe; none leaves the SVM weights as trained.
TEST_F(CrossValidationTest, MergeNormalizationSchemes)
{
    std::map<unsigned int, double> scoreTest, scoreOther;
    std::vector<std::vector<double> > wTest, wOther;
    std::vector<double> featuresTest, featuresOther;
    {
        SCOPED_TRACE("one fold");
        trainAndMerge(1u, MergeNormalization::TEST, scoreTest, wTest,
                      featuresTest);
        trainAndMerge(1u, MergeNormalization::TRAIN, scoreOther, wOther,
                      featuresOther);
        EXPECT_EQ(wTest, wOther);
        ASSERT_EQ(scoreTest.size(), scoreOther.size());
        ASSERT_EQ(900u, scoreTest.size());
        // one PSM per scan, compared in scan order
        auto it = scoreTest.begin();
        auto jt = scoreOther.begin();
        for (; it != scoreTest.end(); ++it, ++jt) {
            EXPECT_EQ(it->second, jt->second);
        }
    }
    for (MergeNormalization scheme :
         {MergeNormalization::TEST, MergeNormalization::TRAIN,
          MergeNormalization::TRAIN_DECOYS, MergeNormalization::NONE}) {
        SCOPED_TRACE(::testing::Message() << "three folds, scheme "
                                          << static_cast<int>(scheme));
        trainAndMerge(3u, scheme, scoreOther, wOther, featuresOther);
        EXPECT_EQ(900u, scoreOther.size());
        EXPECT_EQ(3u, wOther.size());
    }
}

// The pooled direction check only takes the decision. On these data both
// checks keep the trained direction, and the merged scores and weights of
// pooled equal those of bins, with one fold (where pooled is bins) and with
// three. When the check fails, here by a test FDR that no target reaches,
// trainAndMerge finds each test bin scored with its own fold's final
// weights: with --override the trained ones, which were not trained on that
// bin, and otherwise the default vector of a single feature chosen on the
// fold's training set.
TEST_F(CrossValidationTest, DirectionCheckPooledAndFailedChecks)
{
    auto isSingleFeature = [](const std::vector<double>& w) {
        return (w[0] == 0.0) != (w[1] == 0.0);
    };
    for (unsigned int numFolds : {1u, 3u}) {
        for (MergeNormalization scheme :
             {MergeNormalization::TEST, MergeNormalization::TRAIN_DECOYS}) {
            SCOPED_TRACE(::testing::Message() << numFolds << " folds, scheme "
                                              << static_cast<int>(scheme));
            std::map<unsigned int, double> scoreBins, scorePooled;
            std::vector<std::vector<double> > wBins, wPooled;
            std::vector<double> featuresBins, featuresPooled;
            trainAndMerge(numFolds, scheme, scoreBins, wBins, featuresBins);
            DirectionCheckOptions pooled;
            pooled.check = SanityCheck::DirectionCheck::POOLED;
            trainAndMerge(numFolds, scheme, scorePooled, wPooled,
                          featuresPooled, pooled);
            ASSERT_EQ(numFolds, wBins.size());
            for (const std::vector<double>& w : wBins) {
                EXPECT_FALSE(isSingleFeature(w));  // trained, not reset
            }
            EXPECT_EQ(wBins, wPooled);
            ASSERT_EQ(scoreBins.size(), scorePooled.size());
            for (auto it = scoreBins.begin(), jt = scorePooled.begin();
                 it != scoreBins.end(); ++it, ++jt) {
                EXPECT_EQ(it->first, jt->first);
                EXPECT_EQ(it->second, jt->second);
            }
        }
    }
    for (bool overrule : {false, true}) {
        for (bool pooled : {false, true}) {
            for (MergeNormalization scheme :
                 {MergeNormalization::TEST, MergeNormalization::TRAIN_DECOYS,
                  MergeNormalization::NONE}) {
                SCOPED_TRACE(::testing::Message()
                             << "failed check, override " << overrule
                             << ", pooled " << pooled << ", scheme "
                             << static_cast<int>(scheme));
                DirectionCheckOptions failing;
                failing.overrule = overrule;
                failing.check = pooled ? SanityCheck::DirectionCheck::POOLED
                                       : SanityCheck::DirectionCheck::BINS;
                failing.testFdr = 1e-6;
                std::map<unsigned int, double> scoreOf;
                std::vector<std::vector<double> > w;
                std::vector<double> features;
                trainAndMerge(3u, scheme, scoreOf, w, features, failing);
                ASSERT_EQ(3u, w.size());
                for (const std::vector<double>& wFold : w) {
                    EXPECT_EQ(!overrule, isSingleFeature(wFold));
                }
                if (overrule) {
                    // each fold keeps its own trained weights
                    EXPECT_NE(w[0], w[1]);
                    EXPECT_NE(w[0], w[2]);
                }
            }
        }
    }
}

// Whether all folds hold a single-feature vector, as after a reset
static bool allSingleFeature(const std::vector<std::vector<double> >& w) {
    for (const std::vector<double>& wFold : w) {
        if ((wFold[0] == 0.0) == (wFold[1] == 0.0)) {
            return false;
        }
    }
    return true;
}

static void expectSameScores(const std::map<unsigned int, double>& a,
                             const std::map<unsigned int, double>& b) {
    ASSERT_EQ(a.size(), b.size());
    for (auto it = a.begin(), jt = b.begin(); it != a.end(); ++it, ++jt) {
        EXPECT_EQ(it->first, jt->first);
        EXPECT_EQ(it->second, jt->second);
    }
}

// The second chance on small test bins with signal. At --testFDR 0.01 a
// test bin of 300 PSMs would need about 100 targets above its best decoy,
// so the per-bin count fails and bins falls back to single-feature
// vectors, merged with the test scheme, as before the second chance. The
// merged out-of-fold list of 900 PSMs does reach q < 0.01, more often than
// with the initial direction, so second-chance keeps the trained SVMs and
// merges with train-decoys: the same scores and weights as pooled. Another
// scheme set for that case is used instead, and with --override the second
// chance changes nothing. The checks leave the random number generator as
// they found it, also with the pi0 of separate searches, whose estimate
// does draw from it.
TEST_F(CrossValidationTest, SecondChanceKeepsTheSvmsOfSmallBins)
{
    typedef SanityCheck::DirectionCheck Check;
    for (bool separate : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "separate " << separate);
        DirectionCheckOptions options;
        options.testFdr = 0.01;
        options.separate = separate;
        std::map<unsigned int, double> scoreBins, scoreSc, scorePooled;
        std::vector<std::vector<double> > wBins, wSc, wPooled;
        std::vector<double> features;
        DirectionCheckResult bins, sc, pooled;
        options.check = Check::BINS;
        trainAndMerge(3u, MergeNormalization::TEST, scoreBins, wBins,
                      features, options, &bins);
        options.check = Check::SECOND_CHANCE;
        trainAndMerge(3u, MergeNormalization::TEST, scoreSc, wSc, features,
                      options, &sc);
        options.check = Check::POOLED;
        trainAndMerge(3u, MergeNormalization::TRAIN_DECOYS, scorePooled,
                      wPooled, features, options, &pooled);

        EXPECT_TRUE(allSingleFeature(wBins));
        EXPECT_FALSE(bins.secondChanceKept);
        EXPECT_EQ(MergeNormalization::TEST, bins.usedScheme);

        EXPECT_TRUE(sc.secondChanceKept);
        EXPECT_EQ(MergeNormalization::TRAIN_DECOYS, sc.usedScheme);
        ASSERT_EQ(3u, wSc.size());
        for (const std::vector<double>& wFold : wSc) {
            EXPECT_NE(0.0, wFold[0]);
            EXPECT_NE(0.0, wFold[1]);
        }
        EXPECT_EQ(wPooled, wSc);
        expectSameScores(scorePooled, scoreSc);

        for (const DirectionCheckResult* r : {&bins, &sc, &pooled}) {
            EXPECT_TRUE(r->checkKeptRandomState);
            EXPECT_EQ(separate, r->pi0EstimateDraws);
        }
        EXPECT_EQ(pooled.randomStateAfterMerge, sc.randomStateAfterMerge);

        // another scheme for the merge after the second chance
        {
            std::map<unsigned int, double> scoreTrain, scorePooledTrain;
            std::vector<std::vector<double> > wTrain, wPooledTrain;
            DirectionCheckResult train;
            options.check = Check::SECOND_CHANCE;
            options.secondChanceScheme = MergeNormalization::TRAIN;
            trainAndMerge(3u, MergeNormalization::TEST, scoreTrain, wTrain,
                          features, options, &train);
            options.secondChanceScheme = MergeNormalization::TRAIN_DECOYS;
            EXPECT_TRUE(train.secondChanceKept);
            EXPECT_EQ(MergeNormalization::TRAIN, train.usedScheme);
            options.check = Check::POOLED;
            trainAndMerge(3u, MergeNormalization::TRAIN, scorePooledTrain,
                          wPooledTrain, features, options);
            EXPECT_EQ(wPooledTrain, wTrain);
            expectSameScores(scorePooledTrain, scoreTrain);
        }
        // --override: the trained SVMs stay with either check, and the
        // second chance does not change the merge
        {
            std::map<unsigned int, double> scoreOBins, scoreOSc;
            std::vector<std::vector<double> > wOBins, wOSc;
            DirectionCheckResult oBins, oSc;
            options.overrule = true;
            options.check = Check::BINS;
            trainAndMerge(3u, MergeNormalization::TEST, scoreOBins, wOBins,
                          features, options, &oBins);
            options.check = Check::SECOND_CHANCE;
            trainAndMerge(3u, MergeNormalization::TEST, scoreOSc, wOSc,
                          features, options, &oSc);
            options.overrule = false;
            EXPECT_FALSE(oSc.secondChanceKept);
            EXPECT_EQ(MergeNormalization::TEST, oSc.usedScheme);
            EXPECT_EQ(wOBins, wOSc);
            expectSameScores(scoreOBins, scoreOSc);
            EXPECT_EQ(oBins.randomStateAfterMerge, oSc.randomStateAfterMerge);
            EXPECT_FALSE(allSingleFeature(wOSc));
        }
    }
}

// The second chance on inputs without signal, where no target is correct:
// the merged out-of-fold list holds no target with q < --testFDR, the
// trained SVMs are not kept, and second-chance gives exactly the scores,
// weights and final random state of bins, also for separate searches,
// whose pi0 estimates draw from the generator.
TEST_F(CrossValidationTest, SecondChanceRejectsAnInputWithoutSignal)
{
    typedef SanityCheck::DirectionCheck Check;
    for (bool separate : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "separate " << separate);
        DirectionCheckOptions options;
        options.signal = false;
        options.separate = separate;
        std::map<unsigned int, double> scoreBins, scoreSc;
        std::vector<std::vector<double> > wBins, wSc;
        std::vector<double> features;
        DirectionCheckResult bins, sc;
        options.check = Check::BINS;
        trainAndMerge(3u, MergeNormalization::TEST, scoreBins, wBins,
                      features, options, &bins);
        options.check = Check::SECOND_CHANCE;
        trainAndMerge(3u, MergeNormalization::TEST, scoreSc, wSc, features,
                      options, &sc);
        EXPECT_TRUE(allSingleFeature(wBins));
        EXPECT_FALSE(sc.secondChanceKept);
        EXPECT_EQ(MergeNormalization::TEST, sc.usedScheme);
        EXPECT_EQ(wBins, wSc);
        expectSameScores(scoreBins, scoreSc);
        EXPECT_TRUE(sc.checkKeptRandomState);
        EXPECT_EQ(separate, sc.pi0EstimateDraws);
        EXPECT_EQ(bins.randomStateAfterMerge, sc.randomStateAfterMerge);
    }
}

// A perfect separation by the trained SVMs: the targets have 4 as the sum
// of their two features, the decoys less than 2, and neither feature
// separates alone. At --testFDR 0.005 no test bin of about 150 targets
// reaches the test FDR, so the per-bin count fails. The merged out-of-fold
// list puts all 450 targets above all 450 decoys, which the final list's
// postMergeStep rejects as "too good separation" and which ended the run
// after the second chance had kept the SVMs. countMergedTargets now counts
// such a list as unusable (-1): second-chance resets and gives the scores,
// weights and final random state of bins, also for separate searches, and
// pooled resets as well. trainAndMerge completes, so the final merge of
// each check accepts its list.
TEST_F(CrossValidationTest, TooGoodSeparationIsNotKeptByTheMergedCount)
{
    typedef SanityCheck::DirectionCheck Check;
    for (bool separate : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "separate " << separate);
        DirectionCheckOptions options;
        options.perfect = true;
        options.testFdr = 0.005;
        options.separate = separate;
        std::map<unsigned int, double> scoreBins, scoreSc, scorePooled;
        std::vector<std::vector<double> > wBins, wSc, wPooled;
        std::vector<double> features;
        DirectionCheckResult bins, sc, pooled;
        options.check = Check::BINS;
        trainAndMerge(3u, MergeNormalization::TEST, scoreBins, wBins,
                      features, options, &bins);
        options.check = Check::SECOND_CHANCE;
        trainAndMerge(3u, MergeNormalization::TEST, scoreSc, wSc, features,
                      options, &sc);
        options.check = Check::POOLED;
        trainAndMerge(3u, MergeNormalization::TRAIN_DECOYS, scorePooled,
                      wPooled, features, options, &pooled);
        EXPECT_TRUE(allSingleFeature(wBins));
        EXPECT_FALSE(sc.secondChanceKept);
        EXPECT_EQ(MergeNormalization::TEST, sc.usedScheme);
        EXPECT_EQ(wBins, wSc);
        expectSameScores(scoreBins, scoreSc);
        EXPECT_EQ(bins.randomStateAfterMerge, sc.randomStateAfterMerge);
        // the same default vectors, rescaled by another merge scheme
        EXPECT_TRUE(allSingleFeature(wPooled));
        ASSERT_EQ(wBins.size(), wPooled.size());
        for (std::size_t fold = 0; fold < wBins.size(); ++fold) {
            for (std::size_t i = 0; i < 2u; ++i) {
                EXPECT_EQ(wBins[fold][i] > 0.0, wPooled[fold][i] > 0.0);
                EXPECT_EQ(wBins[fold][i] < 0.0, wPooled[fold][i] < 0.0);
            }
        }
        for (const DirectionCheckResult* r : {&bins, &sc, &pooled}) {
            EXPECT_TRUE(r->checkKeptRandomState);
        }
    }
    // the trained SVMs do separate perfectly: with --override they are kept
    // by the bins check, and their merged list is too well separated
    {
        DirectionCheckOptions options;
        options.perfect = true;
        options.testFdr = 0.005;
        options.overrule = true;
        options.check = Check::BINS;
        std::map<unsigned int, double> scoreOf;
        std::vector<std::vector<double> > w;
        std::vector<double> features;
        EXPECT_THROW(trainAndMerge(3u, MergeNormalization::TEST, scoreOf, w,
                                   features, options),
                     MyException);
    }
}

// --static rescales nothing (skipNormalizeScores) and ignores the merge
// scheme: every scheme gives the scores, weights and final random state of
// test, and the merge reports test, which makes -N rescale the scores of all
// PSMs. This holds for separate searches, whose pi0 estimate of each test
// bin under test draws random numbers, and after the second chance kept the
// trained SVMs (--testFDR 0.01), whose merge would otherwise use
// train-decoys or the scheme set for that case.
TEST_F(CrossValidationTest, StaticIgnoresTheMergeScheme)
{
    for (bool separate : {false, true}) {
        for (bool secondChance : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "separate " << separate
                                              << ", second chance "
                                              << secondChance);
            DirectionCheckOptions options;
            options.check = SanityCheck::DirectionCheck::SECOND_CHANCE;
            options.skipNormalizeScores = true;
            options.separate = separate;
            if (secondChance) {
                options.testFdr = 0.01;
            }
            std::map<unsigned int, double> scoreTest;
            std::vector<std::vector<double> > wTest;
            std::vector<double> features;
            DirectionCheckResult test;
            trainAndMerge(3u, MergeNormalization::TEST, scoreTest, wTest,
                          features, options, &test);
            EXPECT_EQ(secondChance, test.secondChanceKept);
            EXPECT_EQ(MergeNormalization::TEST, test.usedScheme);
            EXPECT_FALSE(allSingleFeature(wTest));
            for (MergeNormalization scheme :
                 {MergeNormalization::TRAIN, MergeNormalization::TRAIN_DECOYS,
                  MergeNormalization::NONE}) {
                SCOPED_TRACE(static_cast<int>(scheme));
                options.secondChanceScheme = scheme;
                std::map<unsigned int, double> scoreOf;
                std::vector<std::vector<double> > w;
                DirectionCheckResult other;
                trainAndMerge(3u, scheme, scoreOf, w, features, options,
                              &other);
                EXPECT_EQ(secondChance, other.secondChanceKept);
                EXPECT_EQ(MergeNormalization::TEST, other.usedScheme);
                EXPECT_EQ(wTest, w);
                expectSameScores(scoreTest, scoreOf);
                EXPECT_EQ(test.randomStateAfterMerge,
                          other.randomStateAfterMerge);
            }
        }
    }
}

// Where the per-bin count passes, second-chance is bins: the same scores
// and weights, merged with the test scheme, for one, three and five folds.
// bins itself is the check of r2/fixes, which the failed-check test above
// covers with the default options of trainAndMerge.
TEST_F(CrossValidationTest, SecondChanceIsBinsWhereTheBinsPass)
{
    typedef SanityCheck::DirectionCheck Check;
    for (unsigned int numFolds : {1u, 3u, 5u}) {
        for (bool separate : {false, true}) {
            SCOPED_TRACE(::testing::Message() << numFolds << " folds, "
                                              << "separate " << separate);
            DirectionCheckOptions options;
            options.separate = separate;
            std::map<unsigned int, double> scoreBins, scoreSc;
            std::vector<std::vector<double> > wBins, wSc;
            std::vector<double> features;
            DirectionCheckResult bins, sc;
            options.check = Check::BINS;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreBins,
                          wBins, features, options, &bins);
            options.check = Check::SECOND_CHANCE;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreSc, wSc,
                          features, options, &sc);
            EXPECT_FALSE(allSingleFeature(wBins));
            EXPECT_FALSE(sc.secondChanceKept);
            EXPECT_EQ(MergeNormalization::TEST, sc.usedScheme);
            EXPECT_EQ(wBins, wSc);
            expectSameScores(scoreBins, scoreSc);
            EXPECT_EQ(bins.randomStateAfterMerge, sc.randomStateAfterMerge);
        }
    }
}

// Without a chosen scheme, up to the default three folds merge with test,
// so their outputs do not change; more folds, whose test bins are smaller,
// merge with train-decoys.
TEST_F(CrossValidationTest, DefaultMergeNormalizationOfTheFolds)
{
    for (unsigned int numFolds : {1u, 2u, 3u}) {
        SCOPED_TRACE(numFolds);
        EXPECT_EQ(MergeNormalization::TEST,
                  CrossValidation::defaultMergeNormalization(numFolds));
    }
    for (unsigned int numFolds : {4u, 5u, 10u, 1000u}) {
        SCOPED_TRACE(numFolds);
        EXPECT_EQ(MergeNormalization::TRAIN_DECOYS,
                  CrossValidation::defaultMergeNormalization(numFolds));
    }
}

// After a failed check the folds fall back, unless --override keeps their
// trained vectors: with --init-weights every fold to the weights of the file,
// otherwise each fold to the default vector chosen on its training set, a
// single feature here. trainAndMerge finds each bin scored with its own
// fold's final weights in every case. With -N, the PSMs of the training
// subset keep the score of their test fold only where the folds keep their
// trained vectors: when the check passes and with --override; after a
// reset, also to the weights of the file, all PSMs are scored with the
// averaged weights, and the list of test folds is empty.
TEST_F(CrossValidationTest, FallbackOfAFailedCheckAndTheTestFoldsOfMinusN)
{
    const std::string weightFile =
        ::testing::TempDir() + "perc_fallback_w.txt";
    {
        // normalized weights (skipped), then the raw weights that are read;
        // without a line of feature names, which is then not checked
        std::ofstream out(weightFile);
        out << "9\t9\t9\n"
            << "1\t2\t0\n";
    }
    auto isSingleFeature = [](const std::vector<double>& w) {
        return (w[0] == 0.0) != (w[1] == 0.0);
    };
    for (unsigned int numFolds : {3u, 5u}) {
        for (bool withWeightFile : {false, true}) {
            {
                SCOPED_TRACE(::testing::Message()
                             << numFolds << " folds, weights file "
                             << withWeightFile << ", check passes");
                DirectionCheckOptions passing;
                if (withWeightFile) passing.initWeightFN = weightFile;
                DirectionCheckResult result;
                std::map<unsigned int, double> scoreOf;
                std::vector<std::vector<double> > w;
                std::vector<double> features;
                trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                              features, passing, &result);
                EXPECT_EQ(900u, result.numTestFolds);
            }
            for (bool overrule : {false, true}) {
                SCOPED_TRACE(::testing::Message()
                             << numFolds << " folds, weights file "
                             << withWeightFile << ", failed check, override "
                             << overrule);
                DirectionCheckOptions failing;
                failing.overrule = overrule;
                failing.testFdr = 1e-6;
                if (withWeightFile) failing.initWeightFN = weightFile;
                DirectionCheckResult result;
                std::map<unsigned int, double> scoreOf;
                std::vector<std::vector<double> > w;
                std::vector<double> features;
                trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                              features, failing, &result);
                ASSERT_EQ(numFolds, w.size());
                EXPECT_EQ(overrule ? 900u : 0u, result.numTestFolds);
                if (overrule) {
                    for (unsigned int fold = 1; fold < numFolds; ++fold) {
                        EXPECT_NE(w[0], w[fold]);
                    }
                } else if (withWeightFile) {
                    // the file's vector in every fold, rescaled per bin by
                    // the merge: the same ratio of the feature weights
                    for (const std::vector<double>& wFold : w) {
                        ASSERT_NE(0.0, wFold[1]);
                        EXPECT_NEAR(w[0][0] / w[0][1], wFold[0] / wFold[1],
                                    1e-12);
                        EXPECT_FALSE(isSingleFeature(wFold));
                    }
                } else {
                    for (const std::vector<double>& wFold : w) {
                        EXPECT_TRUE(isSingleFeature(wFold));
                    }
                }
            }
        }
    }
    std::remove(weightFile.c_str());
}

// More folds than spectra leave a fold without PSMs, which
// preIterationSetup reports as an error. With as many folds as spectra of
// equal size, each fold gets one spectrum.
TEST_F(CrossValidationTest, EmptyFoldIsAnError)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet *targets = new DataSet();
    DataSet *decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    std::vector<double> features(2 * 8);
    // 2 spectra with 2 targets and 2 decoys each
    for (int i = 0 ; i < 8 ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = &features[2 * i];
        psm->features[0] = static_cast<double>(i);
        psm->features[1] = static_cast<double>(i % 3);
        psm->scan = 1 + i / 4;
        (i % 2 ? decoys : targets)->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);

    for (unsigned int numFolds : {1u, 2u, 3u}) {
        SCOPED_TRACE(numFolds);
        CrossValidation crossValidation(false, false, 0.01, 0.01, 0.01, 0.0,
                                        0.0, 10, false, 1, false, 1, false,
                                        1.0, numFolds);
        Scores scores(false);
        scores.populateWithPSMs(setHandler);
        EXPECT_EQ(2u, scores.getNumSpectra());
        SanityCheck check;
        if (numFolds <= 2u) {
            EXPECT_NO_THROW(
                crossValidation.preIterationSetup(scores, check, norm.get()));
        } else {
            EXPECT_THROW(
                crossValidation.preIterationSetup(scores, check, norm.get()),
                MyException);
        }
    }
    FeatureNames::resetNumFeatures();
}

namespace {

// line_search as it was before it sorted the breakpoints chunk by chunk:
// one sort of all breakpoints. Also gives the numbers of breakpoints and
// of visited breakpoints.
double fullSortLineSearch(double* w, double* w_bar, double lambda, double* o,
                          double* o_bar, const double* Y, int d, int l,
                          double cpos, double cneg, Delta* deltas,
                          int& numBreakpoints, int& numVisited) {
    int i = 0;
    double omegaL = 0.0;
    double omegaR = 0.0;
    double diff = 0.0;
    for (int i = d; i--;) {
        diff = w_bar[i] - w[i];
        omegaL += w[i] * diff;
        omegaR += w_bar[i] * diff;
    }
    omegaL = lambda * omegaL;
    omegaR = lambda * omegaR;
    double L = omegaL;
    double R = omegaR;
    int ii = 0;
    double d2 = 0.0;
    int p = 0;
    for (i = 0; i < l; i++) {
        diff = Y[i] * (o_bar[i] - o[i]);
        if (Y[i] * o[i] < 1) {
            d2 = ((Y[i] == 1) ? cpos : cneg) * (o_bar[i] - o[i]);
            L += (o[i] - Y[i]) * d2;
            R += (o_bar[i] - Y[i]) * d2;
            if (diff > 0) {
                deltas[p].delta = (1 - Y[i] * o[i]) / diff;
                deltas[p].index = i;
                deltas[p].s = -1;
                p++;
            }
        } else {
            if (diff < 0) {
                deltas[p].delta = (1 - Y[i] * o[i]) / diff;
                deltas[p].index = i;
                deltas[p].s = 1;
                p++;
            }
        }
    }
    std::sort(deltas, deltas + p);
    numBreakpoints = p;
    double delta_prime = 0.0;
    for (i = 0; i < p; i++) {
        numVisited = i;
        delta_prime = L + deltas[i].delta * (R - L);
        if (delta_prime >= 0) {
            break;
        }
        ii = deltas[i].index;
        diff = (deltas[i].s) * ((Y[ii] == 1) ? cpos : cneg) * (o_bar[ii] - o[ii]);
        L += diff * (o[ii] - Y[ii]);
        R += diff * (o_bar[ii] - Y[ii]);
        numVisited = i + 1;
    }
    return (-L / (R - L));
}

}  // namespace

// The chunked sort of the breakpoints visits them in the same order as one
// sort of all of them, so line_search returns the same bits: for random
// problems of up to 40,000 examples, with outputs on a coarse grid so that
// many breakpoints tie, and with searches that stop in the first chunk,
// in a later one, or visit all breakpoints.
TEST(LineSearchTest, ChunkedSortGivesTheSameStepAsAFullSort)
{
#ifdef _OPENMP
    // the parallel-mode sorts of the reference on one thread: hundreds of
    // them are slow on a loaded machine
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    unsigned int seed = 2024u;
    auto uniform = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) + 0.5) / 16777216.0;
    };
    // rounds that visit all breakpoints, only the first chunk, or more
    int numAll = 0, numFirstChunk = 0, numLaterChunks = 0;
    for (int round = 0; round < 200; ++round) {
        const int d = 1 + static_cast<int>(uniform() * 8);
        const int l = round < 20 ? static_cast<int>(uniform() * 50)
                                 : static_cast<int>(uniform() * 40000);
        const bool coarse = round % 3 == 0;
        const double shift = 3.0 * (uniform() - 0.5);
        std::vector<double> w(d), w_bar(d), o(l), o_bar(l), Y(l);
        for (int k = 0; k < d; ++k) {
            w[k] = uniform() - 0.5;
            w_bar[k] = uniform() - 0.5;
        }
        // in every fifth round, the weights shrink and all examples move
        // past the margin, so that the derivative stays negative at every
        // breakpoint and all of them are visited
        const bool visitAll = round % 5 == 1;
        if (visitAll) {
            for (int k = 0; k < d; ++k) w_bar[k] = 0.25 * w[k];
        }
        for (int i = 0; i < l; ++i) {
            Y[i] = uniform() < 0.5 ? 1.0 : -1.0;
            o[i] = 4.0 * (uniform() - 0.5);
            o_bar[i] = o[i] + 2.0 * (uniform() - 0.5) + shift;
            if (visitAll) {
                o[i] = Y[i] * (1.5 * uniform() - 1.0);
                o_bar[i] = Y[i] * (1.0 + uniform());
            }
            if (coarse) {
                o[i] = std::round(o[i] * 4.0) / 4.0;
                o_bar[i] = std::round(o_bar[i] * 4.0) / 4.0;
            }
        }
        const double cpos = 0.1 + uniform() * 10.0;
        const double cneg = 0.1 + uniform() * 10.0;
        std::vector<Delta> deltas(l + 1), referenceDeltas(l + 1);
        const double step = line_search(w.data(), w_bar.data(), 1.0, o.data(),
                                        o_bar.data(), Y.data(), d, l, cpos,
                                        cneg, deltas.data());
        int numBreakpoints = 0, numVisited = 0;
        const double reference = fullSortLineSearch(
            w.data(), w_bar.data(), 1.0, o.data(), o_bar.data(), Y.data(), d,
            l, cpos, cneg, referenceDeltas.data(), numBreakpoints, numVisited);
        if (std::isnan(reference)) {
            EXPECT_TRUE(std::isnan(step)) << "round " << round;
        } else {
            EXPECT_EQ(reference, step) << "round " << round;
        }
        if (numVisited == numBreakpoints) {
            ++numAll;
        } else if (numVisited < std::max(numBreakpoints / 16, 256)) {
            ++numFirstChunk;
        } else {
            ++numLaterChunks;
        }
    }
    EXPECT_GT(numAll, 5);
    EXPECT_GT(numFirstChunk, 5);
    EXPECT_GT(numLaterChunks, 5);
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
}

/* When the trained direction fails the sanity check, the initial vector
 * read with -W is restored, not the default direction (the best single
 * feature), which -W had replaced.
 */
TEST_F(CrossValidationTest, ResetRestoresWeightFileVector)
{
    FeatureNames::resetNumFeatures();
    // the weights file must name the features of the input
    DataSet::resetFeatureNames();
    DataSet::getFeatureNames().insertFeature("f1");
    DataSet::getFeatureNames().insertFeature("f2");
    FeatureNames::setNumFeatures(2);
    const std::string weightFile = ::testing::TempDir() + "perc_reset_w.txt";
    {
        std::ofstream out(weightFile);
        out << "# weights\nf1\tf2\tm0\n"
            << "9\t9\t9\n"           // normalized weights, skipped
            << "0.5\t-2\t0.25\n";    // raw weights, which are read
    }
    // the second feature separates the targets, the first does not
    std::vector<double> features;
    std::vector<PSMDescription> psms(40);
    for (int i = 0; i < 40; ++i) {
        const bool target = i % 2 == 0;
        features.push_back(static_cast<double>(i % 5));
        features.push_back(target ? 10.0 + i : -10.0 - i);
    }
    std::vector<Scores> sets(1, Scores(false));
    for (int i = 0; i < 40; ++i) {
        psms[i].features = &features[2 * i];
        psms[i].scan = static_cast<unsigned int>(i);
        sets[0].addScoreHolder(ScoreHolder(0.0,
            i % 2 == 0 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
    }
    sets[0].recalculateSizes();
    NoNormalizer norm;
    std::vector<double*> none, noneRt;
    norm.setSet(none, noneRt, 2, 0);
    for (bool withWeightFile : {true, false}) {
        SanityCheck check;
        if (withWeightFile) check.setInitWeightFN(weightFile);
        std::vector<std::vector<double> > w(1, std::vector<double>(3, 0.0));
        check.getInitDirection(sets, sets, &norm, w, 0.01, 0.01);
        const std::vector<double> initial = w[0];
        if (withWeightFile) {
            EXPECT_EQ(initial, std::vector<double>({0.5, -2.0, 0.25}));
        }
        // both fallbacks replace the trained bias as well: the default
        // direction sets it to 0, as the initial one
        std::vector<std::vector<double> > trained(1, {7.0, 7.0, 7.0});
        check.resetDirection(sets, trained);
        EXPECT_EQ(trained[0], initial) << "weight file " << withWeightFile;
    }
    std::remove(weightFile.c_str());
    DataSet::resetFeatureNames();
    FeatureNames::resetNumFeatures();
}
