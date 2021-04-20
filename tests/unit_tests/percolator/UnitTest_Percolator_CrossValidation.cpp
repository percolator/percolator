/*
 * Copyright 2020 Brian Raiter <breadbox@muppetlabs.com>
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
 * Unit tests for the CrossValidation class.
 */

#include <gtest/gtest.h>
#include "SetHandler.h"
#include "CrossValidation.h"
#include "Globals.h"

class CrossValidationTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
    CrossValidation *crossValidation;
};

void CrossValidationTest::SetUp()
{
    crossValidation = new CrossValidation(false,  // quickValidation
                                          false,  // reportEachIteration
                                          0.01,   // testFdr
                                          0.01,   // selectionFdr
                                          0.01,   // initialSectionFdr
                                          0.0,    // selectedCpos
                                          0.0,    // selectedCneg
                                          10,     // nIter
                                          true,   // usePi0
                                          1,      // nestedXvalBins
                                          false,  // trainBestPositive
                                          1,      // numThreads
                                          false); // skipNormalizeScores
}

void CrossValidationTest::TearDown()
{
    //delete crossValidation;
}

TEST_F(CrossValidationTest, doStepTest)
{
    int const N = 100;
    int scanNumber = 1;

    Globals::getInstance()->setVerbose(3);

    FeatureMemoryPool featurePool;
    FeatureNames::setNumFeatures(2);

    SetHandler setHandler(0);

    DataSet *targets = new DataSet();
    targets->setLabel(+1);
    for (int ix = 0 ; ix < 2 * N ; ++ix) {
        PSMDescription *psm = new PSMDescription();
        psm->features = new double[2];
        psm->features[0] = static_cast<double>(ix) /
                           static_cast<double>(N) - 1E-6;
        psm->features[1] = static_cast<double>(ix % 2);
        psm->scan = scanNumber++;
        psm->peptide = "T";
        targets->registerPsm(psm);
    }
    setHandler.push_back_dataset(targets);

    DataSet *decoys = new DataSet();
    decoys->setLabel(-1);
    for (int ix = 0 ; ix < N ; ++ix) {
        PSMDescription *psm = new PSMDescription();
        psm->features = new double[2];
        psm->features[0] = static_cast<double>(ix) /
                           static_cast<double>(N);
        psm->features[1] = static_cast<double>(ix % 2);
        psm->scan = scanNumber++;
        psm->peptide = "D";
        decoys->registerPsm(psm);
    }
    setHandler.push_back_dataset(decoys);

    Normalizer *pNorm = NULL;
    setHandler.normalizeFeatures(pNorm);

    SanityCheck *pCheck = new SanityCheck();
    pCheck->checkAndSetDefaultDir();
    pCheck->setConcatenatedSearch(false);

    Scores scores(true);
    scores.populateWithPSMs(setHandler);

    crossValidation->preIterationSetup(scores, pCheck, pNorm, featurePool);

    int found = crossValidation->doStep(true, pNorm, 0.01);
    cerr << "---- got this far: found = " << found << "\n";
}
