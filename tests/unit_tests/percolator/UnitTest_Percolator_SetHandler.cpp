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
#include <sstream>
#include "SetHandler.h"

/* A simple class that tracks global deletions.
 */
class DeletionTracker {
  public:
    DeletionTracker() { }
    virtual ~DeletionTracker() { ++deletionCount; }
    static void reset() { deletionCount = 0; }
    static int deletionCount;
};
int DeletionTracker::deletionCount = 0;

/* A SetHandler subclass with access to some protected methods.
 */
class SetHandlerEx : public SetHandler {
  public:
    SetHandlerEx(int maxPSMs) : SetHandler(maxPSMs) { }
    virtual ~SetHandlerEx() { }
    int getOptionalFieldsEx(const std::string& headerLine, 
                            std::vector<OptionalField>& optionalFields) {
        return this->getOptionalFields(headerLine, optionalFields);
    }
};

/*
 * Unit tests for the SetHandler classes.
 */

class SetHandlerTest : public ::testing::Test {
  protected:
    int testInput(SetHandler *pSH, char const *input);
};

/* For tracking DataSet destructions.
 */
class TrackedDataSet : public DataSet, DeletionTracker { };

// Given a string containing input data, wrap it in a stream object
// and feed it to SetHandler::readTab().
int SetHandlerTest::testInput(SetHandler *handler, char const *input)
{
    std::istringstream str(input);
    SanityCheck *pCheck = NULL;
    int n = handler->readTab(str, pCheck);
    if (n > 0)
        EXPECT_TRUE(pCheck != NULL);
    delete pCheck;
    return n;
}

// Verify DataSets are stored correctly, and are deleted when
// requested.
TEST_F(SetHandlerTest, TestSubsetStoring)
{
    SetHandler sh(0);
    TrackedDataSet *pos = new TrackedDataSet();
    TrackedDataSet *neg = new TrackedDataSet();
    pos->setLabel(+1);
    neg->setLabel(-1);
    DeletionTracker::reset();
    sh.push_back_dataset(pos);
    sh.push_back_dataset(neg);
    EXPECT_EQ(0, DeletionTracker::deletionCount);
    EXPECT_EQ(pos, sh.getSubsetFromLabel(+1));
    EXPECT_EQ(neg, sh.getSubsetFromLabel(-1));

    sh.reset();
    EXPECT_EQ(2, DeletionTracker::deletionCount);
}

// Validate basic parsing of optional fields.
TEST_F(SetHandlerTest, TestGetOptionalFields)
{
    Globals::getInstance()->setVerbose(0);
    SetHandlerEx sh(0);
    vector<OptionalField> fields;
    EXPECT_EQ(0, sh.getOptionalFieldsEx("", fields));
    EXPECT_EQ(0, sh.getOptionalFieldsEx("id\tlabel", fields));
    EXPECT_EQ(1, sh.getOptionalFieldsEx("id\tlabel\tScanNr", fields));
    fields.empty();
    EXPECT_EQ(1, sh.getOptionalFieldsEx("id\tlabel\tfoo\tdm\tbar", fields));
    fields.empty();
    string const realWorldHeader =
        "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tdeltLCn\tdeltCn"
        "\tRefactoredXCorr\tNegLog10PValue\tNegLog10ResEvPValue"
        "\tNegLog10CombinePValue\tPepLen"
        "\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5"
        "\tenzN\tenzC\tenzInt\tlnNumDSP\tdM\tabsdM\tPeptide\tProteins";
    EXPECT_EQ(4, sh.getOptionalFieldsEx(realWorldHeader, fields));
}

// Verify that these input files are rejected.
TEST_F(SetHandlerTest, TestBadReads)
{
    SetHandler sh(0);
    // Empty input.
    EXPECT_THROW(testInput(&sh, ""), MyException);
    // Header line but no data.
    EXPECT_THROW(testInput(&sh, "id\tLabel\tFoo\tPeptide\tProtein\n"),
                 MyException);
    // No feature present.
    EXPECT_THROW(testInput(&sh, "id\tLabel\tPeptide\tProtein\n"
                                "id\t1\t\tp\n"),
                 MyException);
    EXPECT_THROW(testInput(&sh, "id\tLabel\tScanNr\tPeptide\tProtein\n"
                                "id\t1\t0\t\tp\n"),
                 MyException);
    // Missing entry in protein column.
    EXPECT_THROW(testInput(&sh, "id\tLabel\tFeature\tPeptide\tProtein\n"
                                "id\t1\t0\t\t\n"),
                 MyException);
    // Incorrect file type return false instead of throwing
    // (unfortunately it also emits an error message).
    EXPECT_FALSE(testInput(&sh, "<?xml>\n"));
}

// Validate a minimal legal input file.
TEST_F(SetHandlerTest, TestGoodReads)
{
    SetHandler sh(0);
    // Minimum number of columns.
    EXPECT_TRUE(testInput(&sh,
            "id\tLabel\tFeature\tPeptide\tProtein\n"
            "id\t1\t0\t\tP\n"));
    // With a feature column.
    EXPECT_TRUE(testInput(&sh,
            "id\tLabel\tScanNr\tFeature\tPeptide\tProtein\n"
            "id\t1\t0\t0\t\tP\n"));
}

// Validate an input file with a mix of feature columns.
TEST_F(SetHandlerTest, TestReadWithFeatures)
{
    SetHandler sh(0);
    EXPECT_TRUE(testInput(&sh,
            "id\tLabel\tExpMass\tCalcMass\tdM\tabsdM\tPeptide\tProtein\n"
            "id01\t+1\t1324.73\t1324.719\t0.001807\t0.001807\tPEP\tPRO\n"
            "id02\t+1\t3709.60\t3709.646\t-0.01959\t0.019590\tPEP\tPRO\n"
            "id03\t+1\t1998.12\t1998.123\t0.000228\t0.000228\tPEP\tPRO\n"
            "id04\t+1\t3188.54\t3187.612\t-0.00899\t0.008990\tPEP\tPRO\n"
            "id05\t-1\t3838.10\t2837.188\t0.021003\t0.021003\tPEP\tPRO\n"
            "id06\t-1\t2182.15\t2182.175\t-0.02667\t0.026670\tPEP\tPRO\n"));
}
