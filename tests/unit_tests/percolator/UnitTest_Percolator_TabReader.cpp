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
 * Unit tests for the TabReader class.
 */

#include <gtest/gtest.h>
#include <fstream>
#include <random>
#include "DataSet.h"
#include "TmpDir.h"

class TabReaderTest : public ::testing::Test {
  protected:
    void tryValidInt(const std::string& line, int value);
    void tryValidDouble(const std::string& line, double value);
    void tryInvalidInt(const std::string& line);
    void tryInvalidDouble(const std::string& line);
};

// Functions to read a single field.

void TabReaderTest::tryValidInt(const std::string& line, int value)
{
    TabReader reader(line);
    ASSERT_EQ(value, reader.readInt());
    ASSERT_FALSE(reader.error());
}

void TabReaderTest::tryInvalidInt(const std::string& line)
{
    TabReader reader(line);
    reader.readInt();
    ASSERT_TRUE(reader.error());
}

void TabReaderTest::tryValidDouble(const std::string& line, double value)
{
    TabReader reader(line);
    ASSERT_EQ(value, reader.readDouble());
    ASSERT_FALSE(reader.error());
}

void TabReaderTest::tryInvalidDouble(const std::string& line)
{
    TabReader reader(line);
    reader.readDouble();
    ASSERT_TRUE(reader.error());
}

// What is and isn't a valid integer.
TEST_F(TabReaderTest, CheckReadInt)
{
    tryInvalidInt("");
    tryValidInt("-125", -125);
    tryValidInt(" -125", -125);
    tryInvalidInt("- 125");
    tryValidInt("-125 ", -125);
    tryInvalidInt("x125");
    tryInvalidInt("125x");

    tryValidInt("2147483647", 2147483647);
    tryValidInt("-2147483648", -2147483648);
    tryInvalidInt("100000000000000000000000");
    tryInvalidInt("-100000000000000000000000");
}

// What is and isn't a valid real number.
TEST_F(TabReaderTest, CheckReadDouble)
{
    tryInvalidDouble("");
    tryValidDouble("1.5", 1.5);
    tryValidDouble("+1.5\n", +1.5);
    tryValidDouble("-1.5\t", -1.5);
    tryInvalidDouble("x1.5");
    tryInvalidDouble("1:5");
    tryInvalidDouble("15:");
}

// readDouble must give the same value and error state as strtod for every
// input, whichever parsing path it takes.
TEST_F(TabReaderTest, CheckReadDoubleMatchesStrtod)
{
    const char* inputs[] = {
        "0", "-0", "1", "0.1", "-3.14159", "1e10", "1E-10", "2.5e+3",
        "123456789012345678901234567890", "0.1000000000000000055511151231257827",
        "4.9406564584124654e-324", "2.2250738585072014e-308", "1e-400",
        "1.7976931348623157e308", "1e400", "-1e400", "inf", "-inf", "nan",
        "+2.5", " 2.5", "0x1p3", "1.5x", "1.5 ", ".5", "5.", "-", "e5", ""};
    for (const char* input : inputs) {
        std::string line(input);
        char* next = NULL;
        errno = 0;
        double expected = strtod(line.c_str(), &next);
        bool expectError = next == line.c_str() ||
            (*next != '\0' && !isspace(*next)) ||
            ((expected == HUGE_VAL || expected == -HUGE_VAL) && errno == ERANGE);
        TabReader reader(line);
        double value = reader.readDouble();
        EXPECT_EQ(expectError, reader.error() != 0) << "input: " << input;
        if (std::isnan(expected)) {
            EXPECT_TRUE(std::isnan(value)) << "input: " << input;
        } else {
            EXPECT_EQ(0, memcmp(&expected, &value, sizeof(double))) << "input: " << input;
        }
    }
}

// Tests reading a sequence of fields, including skipping fields.
TEST_F(TabReaderTest, CheckSequentialReading)
{
    string test("1\t2\t3\t4\t5");
    TabReader reader(test);
    ASSERT_EQ(1, reader.readInt());
    ASSERT_EQ(2, reader.readInt());
    ASSERT_EQ(3, reader.readInt());
    reader.skip();
    ASSERT_EQ(5, reader.readInt());
    ASSERT_FALSE(reader.error());

    reader.readInt();
    ASSERT_TRUE(reader.error());
}

// The decoy prefix detection as it was before it kept only the label and
// protein columns: a std::string for every column of every row.
static std::string findDecoyPrefixReference(const std::string& fileName,
                                            int proteinIndex, int labelIndex) {
  std::ifstream file(fileName.c_str());
  std::string nextRow;
  getline(file, nextRow);
  std::vector<std::string> proteinNames;
  while (getline(file, nextRow)) {
    TabReader readerRow(nextRow);
    int col = 0;
    bool isDecoy = false;
    while (!readerRow.error()) {
      std::string value = readerRow.readString();
      if (col == labelIndex && value == "-1") isDecoy = true;
      if (col == proteinIndex && isDecoy) proteinNames.push_back(value);
      col++;
    }
  }
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(proteinNames), std::end(proteinNames), rng);
  int numElements = ceil(0.1 * proteinNames.size());
  proteinNames.resize(numElements);
  std::string prefix = TabFileValidator::getLongestCommonPrefix(proteinNames);
  size_t loc = prefix.find("_");
  if (loc != std::string::npos) prefix = prefix.substr(0, loc + 1);
  return prefix;
}

// findDecoyPrefix gives the same prefix as before on generated files, with
// and without common prefixes, with short rows, CRLF line ends and the label
// column after the protein column.
TEST(TabFileValidatorTest, FindDecoyPrefixMatchesReference)
{
  Globals::getInstance()->setVerbose(0);
  TmpDir tmpDir;
  const std::string fileName = tmpDir.getPath() + "/input.pin";
  const char* prefixes[] = {"decoy_", "rev_", "DECOY", "", "XXX_y_", "d"};
  std::mt19937 rng(3);
  for (int trial = 0; trial < 120; ++trial) {
    const int numCols = 3 + static_cast<int>(rng() % 5);
    int labelIndex = static_cast<int>(rng() % numCols);
    int proteinIndex = static_cast<int>(rng() % numCols);
    if (proteinIndex == labelIndex) proteinIndex = (labelIndex + 1) % numCols;
    const char* prefix = prefixes[rng() % 6];
    const size_t numRows = trial < 100 ? rng() % 60 : 20000 + rng() % 20000;
    {
      std::ofstream out(fileName.c_str());
      for (int c = 0; c < numCols; ++c) {
        out << (c ? "\t" : "")
            << (c == labelIndex ? "Label" : c == proteinIndex ? "Proteins" : "f");
      }
      out << "\n";
      for (size_t r = 0; r < numRows; ++r) {
        const bool isDecoy = rng() % 2 == 0;
        const int cols = rng() % 10 == 0 ? static_cast<int>(rng() % numCols) + 1 : numCols;
        for (int c = 0; c < cols; ++c) {
          if (c) out << "\t";
          if (c == labelIndex) {
            out << (isDecoy ? "-1" : rng() % 5 ? "1" : "-1x");
          } else if (c == proteinIndex) {
            if (isDecoy && rng() % 50) out << prefix;
            out << (rng() % 3 ? "prot_" : "p") << rng() % 1000;
          } else {
            out << rng() % 7;
          }
        }
        if (trial % 4 == 1) out << "\r";
        if (r + 1 < numRows || trial % 3) out << "\n";
      }
    }
    EXPECT_EQ(findDecoyPrefixReference(fileName, proteinIndex, labelIndex),
              TabFileValidator::findDecoyPrefix(fileName, proteinIndex, labelIndex))
        << "trial " << trial;
  }
  {
    std::ofstream out(fileName.c_str());
    out << "SpecId\tLabel\tf\tPeptide\tProteins\n";
    for (int r = 0; r < 100; ++r) {
      out << "id" << r << "\t-1\t0.5\tK.PEP.K\tdecoy_prot_" << r << "\tdecoy_x\n";
    }
  }
  EXPECT_EQ("decoy_", TabFileValidator::findDecoyPrefix(fileName, 4, 1));
  EXPECT_EQ("decoy_", TabFileValidator::detectDecoyPrefix(fileName));
  Globals::getInstance()->setVerbose(2);
}
