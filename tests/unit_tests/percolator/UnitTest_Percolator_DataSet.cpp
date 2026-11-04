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
 * Unit tests for the DataSet classes.
 */

#include <gtest/gtest.h>
#include "DataSet.h"

class DataSetTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
    FeatureMemoryPool featurePool;
    PsmStore psmStore;
    std::vector<OptionalField> optionalFields;
    PSMDescription *myPsm;
    std::string decoyPrefix;
};

void DataSetTest::SetUp()
{
    optionalFields.clear();
    featurePool.createPool(1u);
    myPsm = NULL;
    decoyPrefix="dummyVariable_";
}

void DataSetTest::TearDown()
{
    psmStore.clear();
}

// A PSM line needs to at least contain an ID, a Peptide name, and a
// list of Proteins.
TEST_F(DataSetTest, ChecMinimalPsmParsing)
{
    DataSet::readPsm("Id\t-1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ("Id", myPsm->getId());
    ASSERT_EQ("PEPTIDE", myPsm->getFullPeptide());
    ASSERT_EQ(1, myPsm->proteins().size());
    ASSERT_EQ("ProteinList", myPsm->proteins()[0]);
}

// Throw on lines with missing fields.
TEST_F(DataSetTest, CheckInvalidPsmLines)
{
    EXPECT_THROW(DataSet::readPsm("",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
    EXPECT_THROW(DataSet::readPsm("Id",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\t0",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\t0\tPEPTIDE",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
}

TEST_F(DataSetTest, CheckPsmParsingWithFeatures)
{
    optionalFields.push_back(CALCMASS);
    DataSet::readPsm("Id\t1\t1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    // the number of features comes from the first PSM line of the input
    FeatureNames::setNumFeatures(2);
    featurePool.destroyPool();
    featurePool.createPool(2u);
    DataSet::readPsm("Id\t1\t1\t2\t3\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    EXPECT_EQ(3.0, myPsm->features[1]);
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(5);
    featurePool.destroyPool();
    featurePool.createPool(5u);
    DataSet::readPsm("Id\t1\t0.0\t1.1\t2.2\t3.3\t4.4\t5.5\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    EXPECT_EQ(5.5, myPsm->features[4]);
    ASSERT_EQ("PEPTIDE", myPsm->getFullPeptide());
    // a line with one feature column more than the first one
    EXPECT_THROW(DataSet::readPsm(
            "Id\t1\t0.0\t1.1\t2.2\t3.3\t4.4\t5.5\t6.6\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix),
            MyException);
    FeatureNames::resetNumFeatures();
}

// Verify that only integer values are accepted for the label.
TEST_F(DataSetTest, CheckFeatureParsing)
{
    optionalFields.clear();
    DataSet::readPsm("Id\t1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    DataSet::readPsm("Id\t-1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    EXPECT_THROW(DataSet::readPsm("Id\t-1.3\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);

    DataSet::readPsm("Id\t100000000\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    EXPECT_THROW(DataSet::readPsm("Id\t100000000000000000000\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
}

// Verify that features are successfully read.
TEST_F(DataSetTest, CheckOptionalFieldsParsing)
{
    DataSet set;

    optionalFields.push_back(CALCMASS);
    optionalFields.push_back(SCANNR);
    EXPECT_THROW(DataSet::readPsm("Id\t-1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);

    DataSet::readPsm("Id\t-1\t2.5\t42\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ(2.5, myPsm->calcMass);
    ASSERT_EQ(42, myPsm->scan);

    EXPECT_THROW(DataSet::readPsm("Id\t-1\t2:5\t42\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\t-1\t2.5\t4I2\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix), MyException);
}

// The experimental mass enters the comparators of ScoreHolder, so it must be
// a finite number.
TEST_F(DataSetTest, CheckExpMassIsFinite)
{
    optionalFields.push_back(EXPMASS);
    DataSet::readPsm("Id\t1\t1234.5\tPEPTIDE\tProteinList",
            7, optionalFields, true, myPsm, featurePool, psmStore, decoyPrefix);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ(1234.5, myPsm->expMass);
    for (const char* mass : {"nan", "NaN", "inf", "-inf", "Infinity"}) {
        std::string line = std::string("Id\t1\t") + mass + "\tPEPTIDE\tProteinList";
        EXPECT_THROW(DataSet::readPsm(line, 7, optionalFields, true, myPsm,
                featurePool, psmStore, decoyPrefix), MyException) << mass;
    }
}

// A peptide without flanks is used as it is, a number in the peptide column
// (a line with an extra feature column) is an error, and protein inference
// needs the flanks.
TEST_F(DataSetTest, CheckPeptideFormat)
{
    DataSet::readPsm("Id\t1\tK.PEPTIDE.R\tP", 3, optionalFields, true, myPsm,
            featurePool, psmStore, decoyPrefix);
    ASSERT_EQ("PEPTIDE", myPsm->getPeptideSequenceView());
    for (const char* peptide : {"PEPTIDE", "A", "NAN", "INF"}) {
        std::string line = std::string("Id\t1\t") + peptide + "\tP";
        DataSet::readPsm(line, 3, optionalFields, true, myPsm, featurePool,
                psmStore, decoyPrefix);
        EXPECT_EQ(peptide, myPsm->getPeptideSequenceView()) << peptide;
    }
    for (const char* number : {"0.5", "1", "-2e3"}) {
        std::string line = std::string("Id\t1\t") + number + "\tP";
        EXPECT_THROW(DataSet::readPsm(line, 3, optionalFields, true, myPsm,
                featurePool, psmStore, decoyPrefix), MyException) << number;
    }
    ProteinProbEstimator::setCalcProteinLevelProb(true);
    for (const char* peptide : {"PEPTIDE", "K.PEPTIDE", "PEPTIDE.R", "K..R"}) {
        std::string line = std::string("Id\t1\t") + peptide + "\tP";
        EXPECT_THROW(DataSet::readPsm(line, 3, optionalFields, true, myPsm,
                featurePool, psmStore, decoyPrefix), MyException) << peptide;
    }
    DataSet::readPsm("Id\t1\t-.PEPTIDE.R\tP", 3, optionalFields, true, myPsm,
            featurePool, psmStore, decoyPrefix);
    ProteinProbEstimator::setCalcProteinLevelProb(false);
    DataSet::resetWarnings();
}

// The feature pool can be destroyed and created again with another row
// width, and then hands out distinct, writable rows.
TEST(FeatureMemoryPoolTest, ReuseAfterDestroy)
{
    FeatureMemoryPool pool;
    pool.createPool(3);
    double* row = pool.allocate();
    row[0] = row[1] = row[2] = 1.0;
    pool.deallocate(row);
    EXPECT_EQ(row, pool.allocate());  // freed rows are reused
    pool.destroyPool();
    EXPECT_FALSE(pool.isInitialized());

    const size_t width = 7, numRows = 200000;  // more than one block
    pool.createPool(width);
    std::vector<double*> rows;
    for (size_t i = 0; i < numRows; ++i) {
        rows.push_back(pool.allocate());
        for (size_t j = 0; j < width; ++j) rows.back()[j] = static_cast<double>(i);
    }
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < width; ++j) ASSERT_EQ(static_cast<double>(i), rows[i][j]);
    }
}
