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
#include <vector>
#include "SetHandler.h"
#include "CrossValidation.h"
#include "Globals.h"
#include "PseudoRandom.h"

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
    DataSet *targets = new DataSet();
    DataSet *decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    for (int i = 0 ; i < 2 * N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = new double[2];
        psm->features[0] = static_cast<double>(i) / N - 1E-6;
        psm->features[1] = static_cast<double>(i % 2);
        psm->scan = scanNumber++;
        targets->registerPsm(psm);
    }
    for (int i = 0 ; i < N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = new double[2];
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
