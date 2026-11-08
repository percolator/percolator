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
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_map>
#include "SetHandler.h"
#include "ValidateTabFile.h"

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
    SanityCheck check;
    return handler->readTab(str, check);
}

// Verify DataSets are stored correctly, and are deleted when
// requested.
TEST_F(SetHandlerTest, TestSubsetStoring)
{
    SetHandler sh(0);
    TrackedDataSet *pos = new TrackedDataSet();
    TrackedDataSet *neg = new TrackedDataSet();
    pos->setLabel(LabelType::TARGET);
    neg->setLabel(LabelType::DECOY);
    DeletionTracker::reset();
    sh.push_back_dataset(std::unique_ptr<DataSet>(pos));
    sh.push_back_dataset(std::unique_ptr<DataSet>(neg));
    EXPECT_EQ(0, DeletionTracker::deletionCount);
    EXPECT_EQ(pos, sh.getSubsetFromLabel(LabelType::TARGET));
    EXPECT_EQ(neg, sh.getSubsetFromLabel(LabelType::DECOY));

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
    (void)fields.empty();
    EXPECT_EQ(1, sh.getOptionalFieldsEx("id\tlabel\tfoo\tdm\tbar", fields));
    (void)fields.empty();
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

// The Peptide column of the header must follow the features of the first
// PSM line, or the features would be named after other columns.
TEST_F(SetHandlerTest, TestPeptideColumnFollowsFeatures)
{
    SetHandler sh(0);
    // the first PSM has one feature more than the header names
    EXPECT_THROW(testInput(&sh,
            "id\tLabel\tF1\tPeptide\tProtein\n"
            "id\t1\t0\t1\tK.PEP.R\tP\n"),
            MyException);
    sh.reset();
    // the first PSM has a non-numeric feature, so fewer features
    EXPECT_THROW(testInput(&sh,
            "id\tLabel\tF1\tF2\tPeptide\tProtein\n"
            "id\t1\t0\tNA\tK.PEP.R\tP\n"),
            MyException);
    sh.reset();
    EXPECT_TRUE(testInput(&sh,
            "id\tLabel\tF1\tF2\tpeptide\tProtein\n"
            "id\t1\t0\t1\tK.PEP.R\tP\n"
            "id2\t-1\t0\t1\tK.PEP.R\tP\n"));
    sh.reset();
}

// Spectra with more than one PSM of a label are counted, for the warning
// about concatenated searches with several candidates per spectrum, in
// both reading paths; the detection of separate searches is unchanged.
TEST_F(SetHandlerTest, TestCountSpectraWithRepeatedLabels)
{
    const char* input =
        "id\tLabel\tScanNr\tExpMass\tF1\tPeptide\tProtein\n"
        "a\t1\t1\t100\t1\tK.PEP.R\tP\n"
        "b\t-1\t1\t100\t2\tK.PEP.R\tP\n"
        "c\t1\t2\t200\t1\tK.PEP.R\tP\n"
        "d\t1\t2\t200\t3\tK.PEP.R\tP\n"
        "e\t-1\t3\t300\t1\tK.PEP.R\tP\n"
        "f\t-1\t3\t300\t2\tK.PEP.R\tP\n"
        "g\t-1\t3\t300\t2\tK.PEP.R\tP\n";
    for (unsigned int maxPSMs : {0u, 4u}) {
        SetHandler sh(maxPSMs);
        std::istringstream str(input);
        SanityCheck check;
        EXPECT_TRUE(sh.readTab(str, check));
        EXPECT_FALSE(check.concatenatedSearch());
        EXPECT_EQ(7u, sh.getSampledLines());
        EXPECT_EQ(3u, sh.getSampledSpectra());
        EXPECT_EQ(2u, sh.getSampledSpectraWithRepeatedLabel());
        sh.reset();
    }
    SetHandler sh(0);
    std::istringstream str(
        "id\tLabel\tScanNr\tExpMass\tF1\tPeptide\tProtein\n"
        "a\t1\t1\t100\t1\tK.PEP.R\tP\n"
        "b\t-1\t2\t100\t2\tK.PEP.R\tP\n"
        "c\t1\t3\t100\t1\tK.PEP.R\tP\n");
    SanityCheck check;
    EXPECT_TRUE(sh.readTab(str, check));
    EXPECT_TRUE(check.concatenatedSearch());
    EXPECT_EQ(3u, sh.getSampledSpectra());
    EXPECT_EQ(0u, sh.getSampledSpectraWithRepeatedLabel());
    sh.reset();
}

// Line numbers of a concatenated input are mapped back to the input files.
TEST_F(SetHandlerTest, TestDescribeLine)
{
    SetHandler sh(0);
    EXPECT_EQ("line 7", sh.describeLine(7));
    sh.setInputSegments({{2, 2, "a.pin"}, {5, 3, "b.pin"}});
    EXPECT_EQ("line 1", sh.describeLine(1));
    EXPECT_EQ("line 2 of a.pin", sh.describeLine(2));
    EXPECT_EQ("line 4 of a.pin", sh.describeLine(4));
    EXPECT_EQ("line 3 of b.pin", sh.describeLine(5));
    EXPECT_EQ("line 5 of b.pin", sh.describeLine(7));
}

// A file with only a header adds no line to a concatenated input, and the
// segments give the first line of each file with PSMs.
TEST_F(SetHandlerTest, TestConcatenateSkipsFilesWithoutPsms)
{
    namespace fs = std::filesystem;
    fs::path dir = fs::temp_directory_path() / "percolator-unittest-concat";
    fs::create_directories(dir);
    const std::string header = "id\tLabel\tF1\tPeptide\tProteins\n";
    std::vector<std::string> names;
    auto write = [&](const char* name, const std::string& text) {
        std::ofstream((dir / name).string()) << text;
        names.push_back((dir / name).string());
    };
    write("a.pin", header + "a1\t1\t1\tK.PEP.R\tP\na2\t-1\t2\tK.PEP.R\tP\n");
    write("b.pin", header);
    write("c.pin", header + "DefaultDirection\t-\t1\nc1\t1\t3\tK.PEP.R\tP\n");
    std::vector<InputSegment> segments;
    ValidateTabFile validate;
    std::string out = validate.concatenateMultiplePINs(names, dir.string(),
                                                       &segments);
    std::ifstream in(out);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) lines.push_back(line);
    ASSERT_EQ(4u, lines.size());
    EXPECT_EQ("c1\t1\t3\tK.PEP.R\tP\t", lines[3]);
    ASSERT_EQ(2u, segments.size());
    EXPECT_EQ(2u, segments[0].firstLine);
    EXPECT_EQ(2u, segments[0].firstFileLine);
    EXPECT_EQ(names[0], segments[0].fileName);
    EXPECT_EQ(4u, segments[1].firstLine);
    EXPECT_EQ(3u, segments[1].firstFileLine);
    EXPECT_EQ(names[2], segments[1].fileName);
    in.close();
    fs::remove_all(dir);
}

// The search type is detected from scan ids (scan number and experimental
// mass) that occur with both a decoy and another label.
TEST_F(SetHandlerTest, TestSearchTypeDetection)
{
    const char* header = "id\tLabel\tScanNr\tExpMass\tFeature\tPeptide\tProtein\n";
    const char* inputs[][2] = {
        // {input lines, "1" if concatenated}
        {"a\t1\t1\t100.0\t0.5\tK.PEP.K\tP\n"
         "b\t-1\t2\t100.0\t0.5\tK.PEP.K\tP\n", "1"},
        {"a\t1\t1\t100.0\t0.5\tK.PEP.K\tP\n"
         "b\t-1\t1\t100.5\t0.5\tK.PEP.K\tP\n", "1"},
        {"a\t1\t1\t100.0\t0.5\tK.PEP.K\tP\n"
         "b\t-1\t1\t100.0\t0.5\tK.PEP.K\tP\n", "0"},
        {"a\t1\t1\t0.0\t0.5\tK.PEP.K\tP\n"
         "b\t-1\t1\t-0.0\t0.5\tK.PEP.K\tP\n", "0"},
        // a label outside {1,-1} counts as a non-decoy
        {"a\t0\t1\t100.0\t0.5\tK.PEP.K\tP\n"
         "b\t1\t2\t100.0\t0.5\tK.PEP.K\tP\n"
         "c\t-1\t1\t100.0\t0.5\tK.PEP.K\tP\n", "0"}};
    Globals::getInstance()->setVerbose(0);
    for (const auto& input : inputs) {
        for (unsigned int maxPSMs : {0u, 10u}) {  // full read and -N
            SetHandler sh(maxPSMs);
            std::istringstream str(std::string(header) + input[0]);
            SanityCheck check;
            ASSERT_TRUE(sh.readTab(str, check));
            EXPECT_EQ(input[1][0] == '1', check.concatenatedSearch())
                << input[0] << "maxPSMs " << maxPSMs;
        }
    }
    Globals::getInstance()->setVerbose(2);
}

// SeparateSearchDetector against a hash table keyed on (scan, mass), the
// way the search type was detected before.
TEST(SeparateSearchDetectorTest, MatchesHashTable)
{
    struct PairHash {
        size_t operator()(const std::pair<int, double>& p) const {
            return std::hash<int>()(p.first) ^ (std::hash<double>()(p.second) << 1);
        }
    };
    const double masses[] = {0.0, -0.0, 1.5, 2.25, 1e300,
                             std::numeric_limits<double>::quiet_NaN()};
    std::mt19937 rng(7);
    for (int trial = 0; trial < 400; ++trial) {
        // up to 300k scan ids, so that some trials span several checks
        const size_t n = trial < 390 ? rng() % 300 : 100000 + rng() % 200000;
        const int numScans = 1 + static_cast<int>(rng() % (2 * n + 1));
        const bool withDecoys = rng() % 4 != 0;
        std::unordered_map<std::pair<int, double>, bool, PairHash> table;
        bool expected = false;
        SeparateSearchDetector detector;
        bool found = false;
        for (size_t i = 0; i < n; ++i) {
            const int scan = static_cast<int>(rng() % numScans) - numScans / 2;
            const double mass = masses[rng() % 6];
            const bool isDecoy = withDecoys && rng() % 2 == 0;
            auto inserted = table.emplace(std::make_pair(scan, mass), isDecoy);
            if (!inserted.second && inserted.first->second != isDecoy) {
                expected = true;
            }
            const bool foundNow = detector.add(scan, mass, isDecoy);
            EXPECT_TRUE(!found || foundNow);  // stays found
            found = foundNow;
        }
        EXPECT_EQ(expected, detector.foundSeparateSearches()) << "trial " << trial;
    }
}

// A pair of scan ids far apart in the input is found, and a pair within the
// first scan ids is found before the end.
TEST(SeparateSearchDetectorTest, FindsPairsAcrossChecks)
{
    SeparateSearchDetector late;
    late.add(-5, 1234.5, true);
    for (int i = 0; i < 500000; ++i) {
        EXPECT_FALSE(late.add(i, 1234.5, false));
    }
    EXPECT_FALSE(late.add(-5, 1234.5, true));  // the same label again
    EXPECT_FALSE(late.add(-5, 1234.25, false));  // another mass
    late.add(-5, 1234.5, false);
    EXPECT_TRUE(late.foundSeparateSearches());

    SeparateSearchDetector early;
    bool foundEarly = false;
    for (int i = 0; i < 300000 && !foundEarly; ++i) {
        foundEarly = early.add(i / 2, 100.0, i % 2 == 1);
    }
    EXPECT_TRUE(foundEarly);
    EXPECT_TRUE(early.foundSeparateSearches());

    SeparateSearchDetector none;
    for (int i = 0; i < 300000; ++i) {
        none.add(i, 100.0, i % 2 == 1);
    }
    EXPECT_FALSE(none.foundSeparateSearches());
}

// ScanIdPriorityTable against the hash table that the -N path used before:
// the same lookups succeed, with the same priorities, new priorities are
// drawn in the same order, and the byte of each scan id keeps what was
// written to it, here the labels seen so far, also when the table grows.
TEST(ScanIdPriorityTableTest, MatchesHashTable)
{
    struct PairHash {
        size_t operator()(const std::pair<int, double>& p) const {
            return std::hash<int>()(p.first) ^ (std::hash<double>()(p.second) << 1);
        }
    };
    const double masses[] = {0.0, -0.0, 1.5, 2.25, -7.0, 1e300,
                             std::numeric_limits<double>::quiet_NaN()};
    std::mt19937 rng(11);
    for (int trial = 0; trial < 60; ++trial) {
        // the larger trials make the table grow several times
        const size_t n = trial < 50 ? rng() % 500 : 200000 + rng() % 300000;
        const int numScans = 1 + static_cast<int>(rng() % (n + 1));
        std::unordered_map<std::pair<int, double>,
                           std::pair<uint32_t, uint8_t>, PairHash> expected;
        ScanIdPriorityTable table;
        uint32_t counterA = 0, counterB = 0;
        for (size_t i = 0; i < n; ++i) {
            const int scan = static_cast<int>(rng() % numScans) - numScans / 3;
            double mass = masses[rng() % 7];
            if (rng() % 3 == 0) mass += 1.0;
            const uint8_t labelBit = rng() % 2 == 0 ? 1u : 2u;
            auto key = std::make_pair(scan, mass);
            auto it = expected.find(key);
            uint32_t expectedPriority;
            uint8_t expectedByte = 0u;
            if (it != expected.end()) {
                expectedPriority = it->second.first;
                expectedByte = it->second.second;
                it->second.second |= labelBit;
            } else {
                expectedPriority = ++counterA * 2654435761u;
                expected.emplace(key, std::make_pair(expectedPriority, labelBit));
            }
            uint32_t priority = 0;
            bool found = it == expected.end();
            uint8_t& byte = table.findOrInsert(
                scan, mass, [&]() { return ++counterB * 2654435761u; },
                priority, found);
            ASSERT_EQ(it != expected.end(), found) << "trial " << trial << " i " << i;
            ASSERT_EQ(expectedPriority, priority);
            ASSERT_EQ(expectedByte, byte);
            byte |= labelBit;
        }
        EXPECT_EQ(counterA, counterB);
        EXPECT_GE(expected.size(), table.size());  // NaN masses are not stored
    }
}
