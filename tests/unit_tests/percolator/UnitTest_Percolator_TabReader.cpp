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
#include "DataSet.h"

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
