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
    std::vector<OptionalField> optionalFields;
    PSMDescription *myPsm;
};

void DataSetTest::SetUp()
{
    optionalFields.clear();
    featurePool.createPool(1u);
    myPsm = NULL;
}

void DataSetTest::TearDown()
{
    delete myPsm;
}

// A PSM line needs to at least contain an ID, a Peptide name, and a
// list of Proteins.
TEST_F(DataSetTest, ChecMinimalPsmParsing)
{
    DataSet::readPsm("Id\t-1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ("Id", myPsm->getId());
    ASSERT_EQ("PEPTIDE", myPsm->peptide);
    ASSERT_EQ(1, myPsm->proteinIds.size());
    ASSERT_EQ("ProteinList", myPsm->proteinIds[0]);
}

// Throw on lines with missing fields.
TEST_F(DataSetTest, CheckInvalidPsmLines)
{
    EXPECT_THROW(DataSet::readPsm("",
            1, optionalFields, true, myPsm, featurePool), MyException);
    EXPECT_THROW(DataSet::readPsm("Id",
            1, optionalFields, true, myPsm, featurePool), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\t0",
            1, optionalFields, true, myPsm, featurePool), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\t0\tPEPTIDE",
            1, optionalFields, true, myPsm, featurePool), MyException);
    EXPECT_THROW(DataSet::readPsm("Id\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool), MyException);
}

TEST_F(DataSetTest, CheckPsmParsingWithFeatures)
{
    optionalFields.push_back(CALCMASS);
    DataSet::readPsm("Id\t1\t1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
    DataSet::readPsm("Id\t1\t1\t2\t3\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
    DataSet::readPsm("Id\t1\t0.0\t1.1\t2.2\t3.3\t4.4\t5.5\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
}

// Verify that only numerical features are accepted.
TEST_F(DataSetTest, CheckFeatureParsing)
{
    optionalFields.clear();

    DataSet::readPsm("Id\t1\t-1\tPEPTIDE\tProteinList",
        1, optionalFields, true, myPsm, featurePool);
    DataSet::readPsm("Id\t1\t+1.\tPEPTIDE\tProteinList",
        1, optionalFields, true, myPsm, featurePool);
    DataSet::readPsm("Id\t1\t-0.1\tPEPTIDE\tProteinList",
        1, optionalFields, true, myPsm, featurePool);
}


// Verify that features are successfully read.
TEST_F(DataSetTest, CheckOptionalFieldsParsing)
{
    DataSet set;

    optionalFields.push_back(CALCMASS);
    optionalFields.push_back(SCANNR);
    EXPECT_THROW(DataSet::readPsm("Id\t-1\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool), MyException);

    DataSet::readPsm("Id\t-1\t2.5\t42\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ(2.5, myPsm->calcMass);
    ASSERT_EQ(42, myPsm->scan);

    // XXX: If the parser expects an int, it will silently ignore the
    // fractional part of the entry.
    DataSet::readPsm("Id\t-1\t42\t2.5\tPEPTIDE\tProteinList",
            1, optionalFields, true, myPsm, featurePool);
    ASSERT_TRUE(myPsm != NULL);
    ASSERT_EQ(42.0, myPsm->calcMass);
    ASSERT_EQ(2, myPsm->scan);
}
