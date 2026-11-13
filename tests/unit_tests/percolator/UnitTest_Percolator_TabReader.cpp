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
#include <climits>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <locale>
#include <random>
#include <sstream>
#include "DataSet.h"
#include "NumberFormatter.h"
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

// readSimpleDouble, the parser used without std::from_chars, either declines
// a number or gives the same value bits and end as strtod.
TEST_F(TabReaderTest, CheckReadSimpleDoubleMatchesStrtod)
{
  auto check = [](const std::string& text, long& accepted) {
    const char* p = text.c_str();
    double value = 0.0;
    const char* next = NULL;
    if (!TabReader::readSimpleDouble(p, p + text.size(), value, next)) return;
    ++accepted;
    char* strtodNext = NULL;
    const double expected = strtod(p, &strtodNext);
    ASSERT_EQ(0, memcmp(&expected, &value, sizeof(double))) << text;
    ASSERT_EQ(strtodNext, next) << text;
  };
  long accepted = 0;
  const char* edges[] = {
      "0", "-0", "0.0", "-0.000", "1", "-1", "1.", ".5", "-.5", "5.e3",
      "0.1", "3118.1104", "-0.56139", "1e22", "1e-22", "1e23", "1e-23",
      "9007199254740992", "9007199254740993", "1234567890123456789",
      "12345678901234567890", "0.0000000000000000000001",
      "00000000000000000000001.5", "1.0000000000000000000000", "1e",
      "1e+", "1e-", "+1", " 1", "1x", "1 ", "1\t2", "nan", "inf", "-inf",
      "0x10", "1e5x", ".", "-", "-.", "e5", "", "4.9406564584124654e-324",
      "1.7976931348623157e308", "0e99999999", "1e-0", "2.5E+3"};
  for (const char* edge : edges) check(edge, accepted);
  std::mt19937_64 rng(42);
  for (long i = 0; i < 2000000; ++i) {
    std::string text;
    if (rng() & 1) text += '-';
    const int numDigits = 1 + static_cast<int>(rng() % 20);
    const int dot = static_cast<int>(rng() % (numDigits + 2));
    for (int k = 0; k < numDigits; ++k) {
      if (k == dot) text += '.';
      text += static_cast<char>('0' + rng() % 10);
    }
    if (rng() % 4 == 0) {
      text += (rng() & 1) ? 'e' : 'E';
      text += std::to_string(static_cast<int>(rng() % 60) - 30);
    }
    if (rng() % 8 == 0) text += (rng() & 1) ? "\tx" : " ";
    check(text, accepted);
  }
  // most numbers take the fast path where the platform allows it
#if defined(FLT_EVAL_METHOD) && (FLT_EVAL_METHOD == 0 || FLT_EVAL_METHOD == 1)
  EXPECT_GT(accepted, 1000000);
#endif
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

// The decoy prefix detection with a std::string for every column of every
// row: the longest common prefix of every tenth decoy protein, from the
// first one on, up to its first underscore.
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
  std::vector<std::string> sample;
  for (size_t i = 0; i < proteinNames.size(); i += 10) {
    sample.push_back(proteinNames[i]);
  }
  std::string prefix = TabFileValidator::getLongestCommonPrefix(sample);
  size_t loc = prefix.find("_");
  if (loc != std::string::npos) prefix = prefix.substr(0, loc + 1);
  return prefix;
}

// findDecoyPrefix gives the prefix of the reference on generated files, with
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

// The sample is every tenth decoy (the 1st, 11th, 21st, ...), whatever the
// standard library: with decoys of several prefixes the result follows from
// the positions alone. It used to be drawn with std::shuffle and
// std::default_random_engine, which are implementation-defined.
TEST(TabFileValidatorTest, DecoyPrefixSamplesEveryTenthDecoy)
{
  Globals::getInstance()->setVerbose(0);
  TmpDir tmpDir;
  const std::string fileName = tmpDir.getPath() + "/input.pin";
  auto write = [&](int numDecoys, int oddOne) {
    std::ofstream out(fileName.c_str());
    out << "SpecId\tLabel\tf\tPeptide\tProteins\n";
    int decoy = 0;
    for (int r = 0; decoy < numDecoys; ++r) {
      if (r % 3 == 0) {  // targets in between do not count
        out << "t" << r << "\t1\t0.5\tK.PEP.K\tsp_prot" << r << "\n";
        continue;
      }
      // every tenth decoy has the prefix rev_, the others REVERSED_
      const char* prefix = decoy % 10 == 0 ? "rev_" : "REVERSED_";
      if (decoy == oddOne) prefix = "re_";
      out << "d" << r << "\t-1\t0.5\tK.PEP.K\t" << prefix << "prot" << r
          << "\n";
      ++decoy;
    }
  };
  for (int numDecoys : {1, 10, 11, 25, 31, 1000}) {
    write(numDecoys, -1);
    EXPECT_EQ("rev_", TabFileValidator::findDecoyPrefix(fileName, 4, 1))
        << numDecoys << " decoys";
  }
  // the 21st decoy is in the sample of 21 to 30 decoys, and only "re" is
  // common to "re_" and "rev_"; the 2nd decoy is never in it
  write(25, 20);
  EXPECT_EQ("re", TabFileValidator::findDecoyPrefix(fileName, 4, 1));
  write(25, 1);
  EXPECT_EQ("rev_", TabFileValidator::findDecoyPrefix(fileName, 4, 1));
  // no decoys, no prefix
  {
    std::ofstream out(fileName.c_str());
    out << "SpecId\tLabel\tf\tPeptide\tProteins\n"
        << "t1\t1\t0.5\tK.PEP.K\tsp_prot1\n";
  }
  EXPECT_EQ("", TabFileValidator::findDecoyPrefix(fileName, 4, 1));
  Globals::getInstance()->setVerbose(2);
}

// NumberFormatter writes the same characters as operator<< on the stream it
// was made for, on edge values and random bit patterns.
TEST(NumberFormatterTest, MatchesOstream)
{
  std::ostringstream defaultStream;
  EXPECT_EQ(PERCOLATOR_TO_CHARS_DOUBLE != 0,
            NumberFormatter(defaultStream).sameAsStream());
  if (!NumberFormatter(defaultStream).sameAsStream()) {
    GTEST_SKIP() << "std::to_chars for double is not available";
  }
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double edges[] = {
      0.0, -0.0, 1.0, -1.0, 0.5, 2.5, 0.125, 0.1, 1.0 / 3, -2.0 / 3, 1e21,
      1e-5, 1e-4, 9.9999949999999995e-05, 123456, 1234567, 999999.5, 9999995,
      0.99999949999999999, 4.9406564584124654e-324, 2.2250738585072009e-308,
      2.2250738585072014e-308, 1.7976931348623157e308, -1.7976931348623157e308,
      1e100, 1e-100, 3.14159265358979, 1e16, 9007199254740993.0, inf, -inf,
      nan, -nan, std::nan("1")};
  for (int precision = 1; precision <= 17; ++precision) {
    for (double value : edges) {
      std::ostringstream expected;
      expected.precision(precision);
      expected << value;
      const NumberFormatter format(expected);
      ASSERT_TRUE(format.sameAsStream());
      std::string actual;
      format.append(actual, value);
      EXPECT_EQ(expected.str(), actual) << "precision " << precision;
    }
  }
  std::mt19937_64 rng(5);
  for (int i = 0; i < 300000; ++i) {
    uint64_t bits = rng();
    double value;
    memcpy(&value, &bits, sizeof(value));
    if (i % 2) value = static_cast<double>(static_cast<int64_t>(rng() % 2000001) - 1000000) / 100000.0;
    const int precision = i % 5 ? 6 : 1 + static_cast<int>(rng() % 17);
    std::ostringstream expected;
    expected.precision(precision);
    expected << value;
    std::string actual;
    NumberFormatter(expected).append(actual, value);
    ASSERT_EQ(expected.str(), actual) << "bits " << bits << " precision " << precision;
  }
  const long long integers[] = {0, 1, -1, 42, 2147483647, -2147483647 - 1,
                                4294967295LL};
  for (long long value : integers) {
    std::ostringstream expected;
    std::string actual;
    if (value >= INT_MIN && value <= INT_MAX) {
      expected << static_cast<int>(value);
      NumberFormatter::appendInteger(actual, static_cast<int>(value));
      EXPECT_EQ(expected.str(), actual);
    }
    if (value >= 0) {
      expected.str("");
      actual.clear();
      expected << static_cast<unsigned int>(value);
      NumberFormatter::appendInteger(actual, static_cast<unsigned int>(value));
      EXPECT_EQ(expected.str(), actual);
    }
  }
}

// Streams whose output operator<< formats differently are left to operator<<.
TEST(NumberFormatterTest, OtherStreamFormats)
{
  std::vector<std::function<void(std::ostream&)>> settings = {
      [](std::ostream& os) { os << std::fixed; },
      [](std::ostream& os) { os << std::scientific; },
      [](std::ostream& os) { os << std::showpos; },
      [](std::ostream& os) { os << std::showpoint; },
      [](std::ostream& os) { os << std::uppercase; },
      [](std::ostream& os) { os << std::hex; },
      [](std::ostream& os) { os << std::showbase; },
      [](std::ostream& os) { os.width(12); },
      [](std::ostream& os) { os.precision(0); },
      [](std::ostream& os) { os.precision(18); },
      [](std::ostream& os) {
        os.imbue(std::locale(std::locale::classic(), new std::numpunct<char>()));
      }};
  for (size_t i = 0; i < settings.size(); ++i) {
    std::ostringstream os;
    settings[i](os);
    EXPECT_FALSE(NumberFormatter(os).sameAsStream()) << "setting " << i;
  }
  std::ostringstream os;
  os << std::dec << std::left;
  os.precision(10);
  EXPECT_EQ(PERCOLATOR_TO_CHARS_DOUBLE != 0, NumberFormatter(os).sameAsStream());
}
