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
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <unordered_map>
#include "SetHandler.h"
#include "TmpDir.h"
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
    static void rescaleFoldScoresEx(
        Scores& scores,
        const std::vector<std::pair<size_t, unsigned int> >& rows,
        size_t numFolds) {
        rescaleFoldScores(scores, rows, numFolds);
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
    InputSummary summary;
    return handler->readTab(str, summary);
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

// readAndScoreTab scores the PSMs on the lines listed in lineFolds with the
// weights of their fold, and all others with rawWeights (the -N pass).
TEST_F(SetHandlerTest, TestScoreWithFoldWeights)
{
    const char* input =
            "id\tLabel\tScanNr\tf1\tf2\tPeptide\tProtein\n"
            "a\t1\t1\t1.0\t2.0\tK.PEP.K\tP1\n"       // line 2
            "b\t-1\t2\t3.0\t-1.0\tK.PEP.K\tdecoy_P1\n"  // line 3
            "c\t0\t3\t1.0\t1.0\tK.PEP.K\tP1\n"       // line 4, ignored label
            "d\t1\t4\t-2.0\t0.5\tK.PEP.K\tP1\n";     // line 5
    std::vector<double> rawWeights = {1.0, 1.0, 0.5};  // f1, f2, bias
    std::vector<std::vector<double> > foldWeights = {
        {10.0, 0.0, 0.0}, {0.0, 100.0, -1.0}};
    std::vector<std::pair<unsigned int, unsigned int> > lineFolds = {
        {3u, 1u}, {4u, 0u}, {5u, 0u}};
    std::map<std::string, double> expected = {
        {"a", 1.0 + 2.0 + 0.5}, {"b", -100.0 - 1.0}, {"d", -20.0}};
    for (int withFolds = 0; withFolds < 2; ++withFolds) {
        SetHandler sh(0);
        Scores scores(false);
        std::istringstream str(input);
        if (withFolds) {
            ASSERT_TRUE(sh.readAndScoreTab(str, rawWeights, scores,
                                           foldWeights, lineFolds));
        } else {
            ASSERT_TRUE(sh.readAndScoreTab(str, rawWeights, scores));
        }
        ASSERT_EQ(scores.size(), 3u);
        for (const ScoreHolder& h : scores) {
            const std::string id(h.pPSM->getId());
            if (withFolds) {
                EXPECT_DOUBLE_EQ(h.score, expected[id]) << id;
            } else {
                double raw = id == "a" ? 3.5 : (id == "b" ? 2.5 : -1.0);
                EXPECT_DOUBLE_EQ(h.score, raw) << id;
            }
        }
        sh.reset();
    }
}

// The scores of the PSMs scored by a fold get the mean and standard
// deviation of the other PSMs' scores, per fold; the others are unchanged.
TEST_F(SetHandlerTest, TestRescaleFoldScores)
{
    std::vector<PSMDescription> psms(16);
    Scores scores(false);
    const double rest[] = {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -3.0, 0.5};
    std::vector<std::pair<size_t, unsigned int> > rows;
    for (int i = 0; i < 16; ++i) {
        double s;
        if (i < 8) {
            s = rest[i];
        } else if (i < 12) {
            s = 10.0 + 2.0 * (i - 8);  // fold 0: 10, 12, 14, 16
            rows.emplace_back(static_cast<size_t>(i), 0u);
        } else {
            s = -5.0 + 0.5 * (i - 12);  // fold 1: -5, -4.5, -4, -3.5
            rows.emplace_back(static_cast<size_t>(i), 1u);
        }
        scores.addScoreHolder(ScoreHolder(s, LabelType::TARGET, &psms[i]));
    }
    SetHandlerEx::rescaleFoldScoresEx(scores, rows, 2);
    auto moments = [&](int from, int to, double& mean, double& sd) {
        double sum = 0.0, sumSq = 0.0;
        for (int i = from; i < to; ++i) {
            const double s = (scores.begin() + i)->score;
            sum += s;
            sumSq += s * s;
        }
        const double n = to - from;
        mean = sum / n;
        sd = std::sqrt((sumSq - n * mean * mean) / (n - 1.0));
    };
    double restMean, restSd, mean, sd;
    moments(0, 8, restMean, restSd);
    EXPECT_DOUBLE_EQ((scores.begin() + 2)->score, 0.0);  // others unchanged
    EXPECT_NEAR(restMean, 0.0625, 1e-12);
    for (int fold = 0; fold < 2; ++fold) {
        moments(8 + 4 * fold, 12 + 4 * fold, mean, sd);
        EXPECT_NEAR(mean, restMean, 1e-12) << fold;
        EXPECT_NEAR(sd, restSd, 1e-12) << fold;
    }
    // the order within a fold is kept
    EXPECT_LT((scores.begin() + 8)->score, (scores.begin() + 11)->score);
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
        InputSummary summary;
        EXPECT_TRUE(sh.readTab(str, summary));
        EXPECT_FALSE(summary.concatenatedSearch);
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
    InputSummary summary;
    EXPECT_TRUE(sh.readTab(str, summary));
    EXPECT_TRUE(summary.concatenatedSearch);
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
    TmpDir tmpDir;  // removed with its content at the end
    fs::path dir = tmpDir.getPath();
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
    in.close();  // before the directory is removed
}

// Several input files read in turn give the lines of the concatenated file:
// zeros in the columns that a file lacks, no DefaultDirection line but the
// first file's, no line for a file without PSMs. The stream gives the same
// bytes over the borders of its blocks, the segments name the file of each
// line as soon as the line is read, and seekg(0) starts again.
TEST_F(SetHandlerTest, TestSeveralFilesAsOneStream)
{
    namespace fs = std::filesystem;
    TmpDir tmpDir;
    fs::path dir = tmpDir.getPath();
    std::vector<std::string> names;
    auto write = [&](const char* name, const std::string& text) {
        std::ofstream((dir / name).string()) << text;
        names.push_back((dir / name).string());
    };
    Globals::getInstance()->setVerbose(0);
    // a.pin lacks the columns C3 and C4, c.pin lacks C2
    std::string a = "id\tLabel\tF1\tC2\tPeptide\tProteins\n";
    for (int i = 0; i < 5000; ++i) {
        a += "a" + std::to_string(i) + "\t1\t0.5\t1\tK.PEP.R\tP" +
             std::to_string(i) + "\n";
    }
    write("a.pin", a);
    write("b.pin", "id\tLabel\tF1\tC2\tC3\tC4\tPeptide\tProteins\n"
                   "DefaultDirection\t-\t1\t0\t0\t0\n"
                   "b1\t-1\t2\t0\t1\t0\tK.PEP.R\tP\tQ\n");
    write("empty.pin", "id\tLabel\tF1\tC2\tC3\tC4\tPeptide\tProteins\n");
    write("c.pin", "id\tLabel\tF1\tC3\tC4\tPeptide\tProteins\n"
                   "c1\t1\t3\t0\t1\tK.PEP.R\tP\n");
    ::testing::internal::CaptureStderr();
    ConcatenatedPinStream stream(names);
    EXPECT_EQ("Warning: the input file " + names[2] + " contains no PSMs.\n",
              ::testing::internal::GetCapturedStderr());
    EXPECT_TRUE(stream.segments().empty());

    std::vector<std::string> lines;
    std::string line;
    ASSERT_TRUE(static_cast<bool>(std::getline(stream, line)));
    EXPECT_EQ("id\tLabel\tF1\tC2\tC3\tC4\tPeptide\tProteins\t", line);
    lines.push_back(line);
    ASSERT_TRUE(static_cast<bool>(std::getline(stream, line)));
    EXPECT_EQ("a0\t1\t0.5\t1\t0\t0\tK.PEP.R\tP0\t", line);
    lines.push_back(line);
    // the file of a line is known when the line has been read
    ASSERT_EQ(1u, stream.segments().size());
    EXPECT_EQ(2u, stream.segments()[0].firstLine);
    while (std::getline(stream, line)) lines.push_back(line);
    ASSERT_EQ(5003u, lines.size());
    EXPECT_EQ("a4999\t1\t0.5\t1\t0\t0\tK.PEP.R\tP4999\t", lines[5000]);
    EXPECT_EQ("b1\t-1\t2\t0\t1\t0\tK.PEP.R\tP\tQ\t", lines[5001]);
    EXPECT_EQ("c1\t1\t3\t0\t0\t1\tK.PEP.R\tP\t", lines[5002]);
    const std::vector<InputSegment> segments = stream.segments();
    ASSERT_EQ(3u, segments.size());
    EXPECT_EQ(5002u, segments[1].firstLine);
    EXPECT_EQ(3u, segments[1].firstFileLine);  // after its DefaultDirection
    EXPECT_EQ(names[1], segments[1].fileName);
    EXPECT_EQ(5003u, segments[2].firstLine);
    EXPECT_EQ(2u, segments[2].firstFileLine);
    EXPECT_EQ(names[3], segments[2].fileName);

    // the file that concatenateMultiplePINs writes holds the same lines
    std::vector<InputSegment> fileSegments;
    ValidateTabFile validate;
    Globals::getInstance()->setVerbose(0);
    ::testing::internal::CaptureStderr();
    std::string out = validate.concatenateMultiplePINs(names, dir.string(),
                                                       &fileSegments);
    ::testing::internal::GetCapturedStderr();
    {
        std::ifstream in(out);
        std::size_t ix = 0;
        for (; std::getline(in, line); ++ix) {
            ASSERT_LT(ix, lines.size());
            ASSERT_EQ(lines[ix], line) << "line " << ix + 1;
        }
        EXPECT_EQ(lines.size(), ix);
    }
    ASSERT_EQ(segments.size(), fileSegments.size());
    for (std::size_t ix = 0; ix < segments.size(); ++ix) {
        EXPECT_EQ(segments[ix].firstLine, fileSegments[ix].firstLine);
        EXPECT_EQ(segments[ix].firstFileLine, fileSegments[ix].firstFileLine);
        EXPECT_EQ(segments[ix].fileName, fileSegments[ix].fileName);
    }

    // again from the start, after the end and from the middle of a file
    for (int pass = 0; pass < 2; ++pass) {
        stream.clear();
        stream.seekg(0, std::ios::beg);
        ASSERT_TRUE(static_cast<bool>(stream));
        const std::size_t numLines = pass == 0 ? 2500u : lines.size();
        for (std::size_t ix = 0; ix < numLines; ++ix) {
            ASSERT_TRUE(static_cast<bool>(std::getline(stream, line)));
            ASSERT_EQ(lines[ix], line) << "pass " << pass << ", line " << ix + 1;
        }
        EXPECT_EQ(3u, stream.segments().size());
    }
    EXPECT_FALSE(static_cast<bool>(std::getline(stream, line)));
    // no other position can be sought
    stream.clear();
    stream.seekg(10, std::ios::beg);
    EXPECT_TRUE(stream.fail());
}

// The reader of all PSMs and the subset reader of -N give the same
// messages, in the same order, and the same error for malformed lines: the
// error of the lowest failing line, after the messages of the lines before
// it. A reader that works ahead of the line it reports (number conversion
// on several threads was tried and measured slower) has to keep this.
TEST_F(SetHandlerTest, TestBothReadersReportMalformedLinesAlike)
{
    const std::string header = "id\tLabel\tF1\tF2\tPeptide\tProteins\n";
    const int numLines = 9000;
    // (line of the input, what replaces its PSM line)
    typedef std::vector<std::pair<int, std::string> > Edits;
    const std::string strange = "x\t1\t0.5\tnan\tK.PEP.R\tP";
    const std::string text = "x\t1\tabc\t1\tK.PEP.R\tP";
    const std::string ignored = "x\t0\t0.5\t1\tK.PEP.R\tP";
    const std::string unflanked = "x\t1\t0.5\t1\tPEPTIDE\tP";
    const std::string noPeptide = "x\t1\t0.5\t1";
    const std::string blank = "x\t1\t0.5 \t1\tK.PEP.R\tP";  // still read
    const std::string shifted = "x\t1\t0.5 2\t1\tK.PEP.R\tP";
    const std::string strangeAndUnflanked = "x\t1\tinf\t1\tPEPTIDE\tP";
    const std::vector<Edits> cases = {
        {},
        {{10, strange}},
        {{10, ignored}, {20, strange}, {30, ignored}},
        {{10, strange}, {20, ignored}, {30, text}, {40, ignored}},
        {{2, unflanked}, {3, strange}},
        {{2, strange}, {3, unflanked}},
        {{5000, noPeptide}, {10, strange}},
        {{10, noPeptide}, {9, strange}, {11, strange}},
        {{4097, strange}, {4098, ignored}, {8193, text}},
        {{4096, text}, {4097, strange}},
        {{100, blank}, {200, shifted}, {150, ignored}},
        {{7, strangeAndUnflanked}, {8, ignored}},
        {{numLines + 1, strange}, {numLines, ignored}},
        {{2, text}}};
    const int verbosity = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(1);
    int numErrors = 0, numMessages = 0;
    for (std::size_t c = 0; c < cases.size(); ++c) {
        std::vector<std::string> lines(static_cast<std::size_t>(numLines));
        for (int i = 0; i < numLines; ++i) {
            lines[static_cast<std::size_t>(i)] =
                "p" + std::to_string(i) + "\t" + (i % 2 ? "1" : "-1") + "\t" +
                std::to_string(i % 7) + ".25\t1e-3\tK.PEP.R\tP";
        }
        for (const auto& edit : cases[c]) {
            lines[static_cast<std::size_t>(edit.first - 2)] = edit.second;
        }
        std::string input = header;
        for (const std::string& line : lines) input += line + "\n";
        for (bool noTerminate : {false, true}) {
            Globals::getInstance()->setNoTerminate(noTerminate);
            std::string results[2];
            for (int lineByLine = 0; lineByLine < 2; ++lineByLine) {
                DataSet::resetWarnings();
                FeatureNames::resetNumFeatures();
                DataSet::resetFeatureNames();
                SetHandler sh(lineByLine ? 1000000u : 0u);
                std::istringstream str(input);
                InputSummary summary;
                std::string error;
                ::testing::internal::CaptureStderr();
                try {
                    sh.readTab(str, summary);
                } catch (const MyException& e) {
                    error = e.what();
                }
                results[lineByLine] =
                    ::testing::internal::GetCapturedStderr() + "error: " +
                    error;
                if (!lineByLine) {
                    numErrors += error.empty() ? 0 : 1;
                    numMessages += results[0].size() > 7u ? 1 : 0;
                }
            }
            EXPECT_EQ(results[1], results[0])
                << "case " << c << ", --no-terminate " << noTerminate;
        }
    }
    EXPECT_GT(numErrors, 10);
    EXPECT_GT(numMessages, 20);
    Globals::getInstance()->setNoTerminate(false);
    Globals::getInstance()->setVerbose(verbosity);
    DataSet::resetWarnings();
    DataSet::resetFeatureNames();
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
            InputSummary summary;
            ASSERT_TRUE(sh.readTab(str, summary));
            EXPECT_EQ(input[1][0] == '1', summary.concatenatedSearch)
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

/* Gives the tests access to the two ways of reading a PSM line's header.
 */
class SetHandlerHeaderEx : public SetHandler {
  public:
    SetHandlerHeaderEx() : SetHandler(0) { }
    ScanId getScanIdEx(const std::string& line, int& label,
                       std::vector<OptionalField>& fields, unsigned int lineNr) {
        return getScanId(line, label, fields, lineNr);
    }
    static void incrementLineNrEx(unsigned int& lineNr) {
        incrementLineNr(lineNr);
    }
    bool readScanIdAndLabelEx(const std::string& line, unsigned int lineNr,
                              std::vector<OptionalField>& fields, TabReader& reader,
                              PsmHeader& header, ScanId& scanId, int& label) {
        return readScanIdAndLabel(line, lineNr, fields, reader, header, scanId, label);
    }
};

static std::string bitsOf(double value) {
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return std::to_string(bits);
}

static std::string describePsm(PSMDescription* psm, LabelType label) {
    std::ostringstream out;
    out << " label " << static_cast<int>(label) << " id " << psm->getId()
        << " scan " << psm->scan << " exp " << bitsOf(psm->expMass)
        << " calc " << bitsOf(psm->calcMass)
        << " rt " << bitsOf(psm->getRetentionTime())
        << " file " << psm->getSpectrumFileName() << " line " << psm->getLineNr()
        << " pep " << psm->getFullPeptide() << " prot ";
    psm->printProteins(out);
    for (size_t i = 0; i < FeatureNames::getNumFeatures(); ++i) {
        out << " f " << bitsOf(psm->features[i]);
    }
    return out.str();
}

// Reading a PSM line with the fast path of DataSet::readPsmHeader and
// readPsmBody gives the same scan id, label, PSM, messages and exceptions as
// getScanId followed by DataSet::readPsm, on valid and malformed lines.
TEST(PsmHeaderTest, FastPathMatchesGetScanIdAndReadPsm)
{
    const int verbosity = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(1);
    FeatureNames::setNumFeatures(2);
    FeatureMemoryPool featurePool;
    featurePool.createPool(2);
    PsmStore psmStore;
    const std::string decoyPrefix = "decoy_";
    SetHandlerHeaderEx sh;
    const char* numbers[] = {"123", "1.5", "-0.0", "1e3", "0", "-7", "x", "",
        "1.5 ", " 2", "1.5x", "nan", "inf", "1e400", "2147483648", "+4", "1."};
    const char* labels[] = {"1", "-1", "0", "+1", "x", "1 ", " -1", "", "2"};
    const char* texts[] = {"a.raw", "b.raw", "", "c d"};
    const OptionalField allFields[] = {SCANNR, EXPMASS, CALCMASS, RETTIME,
                                       FILENAME, DELTAMASS};
    std::mt19937 rng(17);
    int numSimple = 0;
    // lines that the fast path refuses and that are read all the same, with
    // an example of each kind
    int numRefusedButRead = 0;
    std::map<std::string, std::string> refusedButRead;
    for (int trial = 0; trial < 30000; ++trial) {
        std::vector<OptionalField> fields;
        const int numFields = static_cast<int>(rng() % 5);
        for (int k = 0; k < numFields; ++k) {
            fields.push_back(allFields[rng() % (trial % 50 == 0 ? 6 : 5)]);
        }
        auto pick = [&](const char* const* values, size_t n, size_t good) {
            // mostly well-formed values
            return std::string(values[rng() % 4 ? rng() % good : rng() % n]);
        };
        std::string line = (rng() % 20 ? "id" + std::to_string(trial) : "");
        line += "\t" + pick(labels, 9, 3);
        for (OptionalField field : fields) {
            line += "\t";
            line += field == FILENAME ? pick(texts, 4, 2) : pick(numbers, 17, 6);
        }
        for (int k = 0; k < 2; ++k) line += "\t" + pick(numbers, 17, 6);
        line += rng() % 10 ? "\tK.PEPTIDE.R" : "\t";
        line += rng() % 3 ? "\tP1" : "\tP1\tdecoy_P2";
        if (rng() % 8 == 0) line = line.substr(0, rng() % (line.size() + 1));
        if (rng() % 10 == 0) line.insert(rng() % (line.size() + 1), rng() % 2 ? "\t" : " ");
        const unsigned int lineNr = 3 + static_cast<unsigned int>(trial);
        Globals::getInstance()->setNoTerminate(trial % 3 == 0);

        std::string results[2];
        for (int fast = 0; fast < 2; ++fast) {
            // both ways print the warnings that are printed once per run
            DataSet::resetWarnings();
            std::ostringstream messages;
            std::streambuf* cerrBuffer = std::cerr.rdbuf(messages.rdbuf());
            std::ostringstream result;
            bool simple = false, bodyRefused = false, read = false;
            try {
                int label = 0;
                ScanId scanId;
                TabReader reader(line);
                PsmHeader header;
                if (fast) {
                    simple = sh.readScanIdAndLabelEx(line, lineNr, fields, reader,
                                                     header, scanId, label);
                    numSimple += simple ? 1 : 0;
                } else {
                    scanId = sh.getScanIdEx(line, label, fields, lineNr);
                }
                result << "scan " << scanId.first << " " << bitsOf(scanId.second)
                       << " label " << label;
                if (label == 1 || label == -1) {
                    PSMDescription* psm = NULL;
                    LabelType psmLabel = LabelType::UNDEFINED;
                    if (!simple || !DataSet::readPsmBody(reader, header, true, psm,
                            psmLabel, featurePool, psmStore, decoyPrefix)) {
                        bodyRefused = simple;
                        psmLabel = DataSet::readPsm(line, lineNr, fields, true, psm,
                            featurePool, psmStore, decoyPrefix);
                    }
                    result << describePsm(psm, psmLabel);
                    read = true;
                }
            } catch (const MyException& e) {
                result << "exception " << e.what();
            }
            std::cerr.rdbuf(cerrBuffer);
            results[fast] = result.str() + " messages " + messages.str();
            if (fast && read && (!simple || bodyRefused)) {
                // The kinds listed at SetHandler::PsmLineReader: the header
                // is refused for a non-finite experimental mass, which is
                // only read under --no-terminate, and for a repeated
                // column; the body for a calculated mass or retention time
                // that is not one number ending at its tab.
                ++numRefusedButRead;
                std::string kind = "body: mass or retention time column";
                if (!simple) {
                    kind = messages.str().find("experimental mass") !=
                                   std::string::npos
                               ? "header: non-finite experimental mass"
                               : "header: repeated column";
                    if (kind == "header: repeated column") {
                        std::vector<OptionalField> sorted(fields);
                        std::sort(sorted.begin(), sorted.end());
                        EXPECT_TRUE(std::adjacent_find(sorted.begin(),
                                                       sorted.end()) !=
                                    sorted.end())
                            << "refused without a repeated column: " << line;
                    } else {
                        EXPECT_TRUE(NO_TERMINATE) << line;
                    }
                } else {
                    EXPECT_TRUE(std::find(fields.begin(), fields.end(),
                                          CALCMASS) != fields.end() ||
                                std::find(fields.begin(), fields.end(),
                                          RETTIME) != fields.end())
                        << line;
                }
                refusedButRead.emplace(kind, line);
            }
        }
        ASSERT_EQ(results[0], results[1]) << "line " << line;
    }
    EXPECT_GT(numSimple, 5000);  // the fast path is exercised
    // every kind of refused line that is read occurs
    EXPECT_EQ(3u, refusedButRead.size());
    EXPECT_GT(numRefusedButRead, 100);
    if (::testing::GTEST_FLAG(also_run_disabled_tests)) {  // to list them
        std::cout << numRefusedButRead << " refused lines were read\n";
        for (const auto& example : refusedButRead) {
            std::cout << example.first << ": [" << example.second << "]\n";
        }
    }
    Globals::getInstance()->setNoTerminate(false);
    Globals::getInstance()->setVerbose(verbosity);
    FeatureNames::resetNumFeatures();
    PSMDescription::resetSpectrumFileNames();
}

// A scan number above 2^31 - 1 is refused; the error said "Check if scan
// number is an integer", which it is.
TEST(PsmHeaderTest, ScanNumberErrorNamesTheRange)
{
    SetHandlerHeaderEx sh;
    std::vector<OptionalField> fields = {SCANNR};
    FeatureNames::resetNumFeatures();
    const std::string range = "between -2147483648 and 2147483647.";
    for (const char* scan : {"2147483648", "-2147483649", "4294967296", "1.5", "x"}) {
        const std::string line = std::string("id\t1\t") + scan + "\tK.PEPTIDEK.A\tprot";
        int label = 0;
        try {
            sh.getScanIdEx(line, label, fields, 7u);
            ADD_FAILURE() << "no exception for " << scan;
        } catch (const MyException& e) {
            const std::string message = e.what();
            EXPECT_NE(std::string::npos, message.find("scan number on line 7"))
                << message;
            EXPECT_NE(std::string::npos, message.find(range)) << message;
        }
        FeatureMemoryPool pool;
        PsmStore store;
        PSMDescription* psm = NULL;
        try {
            DataSet::readPsm(line, 7u, fields, true, psm, pool, store, "decoy_");
            ADD_FAILURE() << "no exception for " << scan;
        } catch (const MyException& e) {
            const std::string message = e.what();
            EXPECT_NE(std::string::npos, message.find("scan number of PSM id"))
                << message;
            EXPECT_NE(std::string::npos, message.find(range)) << message;
        }
    }
    // the ends of the range are read
    for (const char* scan : {"2147483647", "-2147483648", "0"}) {
        const std::string line = std::string("id\t1\t") + scan + "\tK.PEPTIDEK.A\tprot";
        int label = 0;
        ScanId scanId;
        EXPECT_NO_THROW(scanId = sh.getScanIdEx(line, label, fields, 7u)) << scan;
        EXPECT_EQ(std::string(scan), std::to_string(scanId.first));
    }
    FeatureNames::resetNumFeatures();
}

// Line numbers and PSM counts are 32-bit integers, also signed ones: line
// 2^31 used to become scan number -2147483648, and counts wrapped.
TEST(PsmHeaderTest, MoreThanIntMaxLinesIsAnError)
{
    unsigned int lineNr = 2u;
    SetHandlerHeaderEx::incrementLineNrEx(lineNr);
    EXPECT_EQ(3u, lineNr);

    lineNr = 2147483646u;
    EXPECT_NO_THROW(SetHandlerHeaderEx::incrementLineNrEx(lineNr));
    EXPECT_EQ(2147483647u, lineNr);  // the last line that can be read
    try {
        SetHandlerHeaderEx::incrementLineNrEx(lineNr);
        FAIL() << "no exception";
    } catch (const MyException& e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find(
            "cannot process more than 2147483647 PSMs in one run")) << e.what();
    }
    EXPECT_EQ(2147483647u, lineNr);
}
