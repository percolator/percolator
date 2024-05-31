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
 * Unit tests for the Option class.
 */

#include <gtest/gtest.h>
#include "Option.h"

class OptionTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    CommandLineParser parser;
    void parseArgv(char const *argv[]);
};

void OptionTest::SetUp()
{
    parser.defineOption("a", "aShortOption",
                        "This is a short option.", "", TRUE_IF_SET);
    parser.defineOption(Option::NO_SHORT_OPT, "longOption",
                        "This is a long option.", "", FALSE_IF_SET);
    parser.defineOption(Option::EXPERIMENTAL_FEATURE, "experimental",
                        "An experimental option.", "", TRUE_IF_SET);
    parser.defineOption("i", "integerOption", "This takes an integer.",
                        "value", "42");
    parser.defineOption("n", "numberOption", "This takes a real number.",
                        "value", "3.1415926");
    parser.defineOption("s", "stringOptions", "This takes anything.",
                        "value", "foo");
}

// Compute argc and call parser.parseArgs().
void OptionTest::parseArgv(char const *argv[])
{
    int argc;
    for (argc = 0 ; argv[argc] != NULL ; ++argc) ;
    parser.parseArgs(argc, const_cast<char**>(argv));
}

TEST_F(OptionTest, CheckBasicOptionParsing)
{
    char const *argv[] = {
        "prog", "-a", "--longOption", "-i", "15", "-n", "1.5", NULL
    };

    parseArgv(argv);

    EXPECT_TRUE(parser.isOptionSet("aShortOption"));
    EXPECT_TRUE(parser.isOptionSet("longOption"));
    EXPECT_TRUE(parser.isOptionSet("integerOption"));
    EXPECT_TRUE(parser.isOptionSet("numberOption"));
    EXPECT_FALSE(parser.isOptionSet("stringOption"));
    EXPECT_FALSE(parser.isOptionSet("experimental"));

    EXPECT_EQ("1", parser.options["aShortOption"]);
    EXPECT_EQ("0", parser.options["longOption"]);
    EXPECT_EQ(15, parser.getInt("integerOption", 0, 100));
    EXPECT_EQ(1.5, parser.getDouble("numberOption", 0.0, 100.0));

    // Default values aren't set by the parser.
    EXPECT_EQ("", parser.options["stringOptions"]);
    // You can't look up options by their short names.
    EXPECT_FALSE(parser.isOptionSet("a"));
}

TEST_F(OptionTest, CheckThrowOnInvalidOption)
{
    char const *argv[] = {
        "prog", "--unknownOption", "foo", NULL
    };

    EXPECT_THROW(parseArgv(argv), MyException);
}

TEST_F(OptionTest, CheckThrowOnRangeExceeded)
{
    char const *argv[] = {
        "prog", "--integerOption", "42", NULL
    };

    parseArgv(argv);

    EXPECT_THROW(parser.getInt("integerOption", 0, 10), MyException);
    EXPECT_EQ(42, parser.getInt("integerOption", 10, 50));
    EXPECT_THROW(parser.getInt("integerOption", 50, 60), MyException);
}

TEST_F(OptionTest, CheckSeparationOfOptionsAndArguments)
{
    char const *argv[] = {
        "prog", "a1", "-a", "a2", "-i", "0", "a3", "-s", "0", "a4", NULL
    };

    parseArgv(argv);

    EXPECT_EQ(4, parser.arguments.size());
    EXPECT_EQ("a1", parser.arguments[0]);
    EXPECT_EQ("a2", parser.arguments[1]);
    EXPECT_EQ("a3", parser.arguments[2]);
    EXPECT_EQ("a4", parser.arguments[3]);
}
