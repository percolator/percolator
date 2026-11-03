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
#include <cmath>
#include <cstdarg>
#include "SetHandler.h"
#include "DataSet.h"
#include "PosteriorEstimator.h"
#include "Scores.h"

// Some strings in alphabetical order.
static std::string const psmNames[] = { "ABC", "DEF", "GHI", "JKL", "MNO" };

class ScoreHolderTest : public ::testing::Test {
  protected:
    bool checkOrder(std::vector<ScoreHolder> const *scores, ...);
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
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
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
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
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
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
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
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
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
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
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
        psm = psmStore.add(PSMDescription(psmNames[i]));
        psm->scan = 100 + i;
        set1->registerPsm(psm);
        psm = psmStore.add(PSMDescription(psmNames[i]));
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
        std::vector<std::pair<double, bool> > combined;
        // countTargetsBelowFdr needs the q values of the decoys as well;
        // calcScoresAndQvals sets this only after the first call below
        PosteriorEstimator::setNegative(true);
        for (double fdr : {0.01, 0.05, 0.2}) {
            for (bool skipDecoysPlusOne : {false, true}) {
                int count = scores.countTargetsBelowFdr(
                    rows, w, fdr, skipDecoysPlusOne, combined);
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
