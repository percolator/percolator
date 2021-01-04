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
 * Unit tests for the PSMDescription and PSMDescriptionDOC classes.
 */

#include <gtest/gtest.h>
#include "PSMDescription.h"

TEST(PSMDescriptionTest, CheckRemovePTMs)
{
    ASSERT_EQ("K.DAQLLVAER.M", PSMDescription::removePTMs("K.DAQLLVAER.M"));
    ASSERT_EQ("A.CDEF.H", PSMDescription::removePTMs("ABCD[0]EFGH"));
    ASSERT_EQ("A.BCDEFGHJ.K", PSMDescription::removePTMs("A!B[123]CDE[456]FGH[789]J?K"));

    EXPECT_THROW(PSMDescription::removePTMs(""), MyException);
    EXPECT_THROW(PSMDescription::removePTMs("1"), MyException);
    EXPECT_THROW(PSMDescription::removePTMs("12"), MyException);
    EXPECT_THROW(PSMDescription::removePTMs("123"), MyException);
    ASSERT_EQ("1..4", PSMDescription::removePTMs("1234"));
}

TEST(PSMDescriptionTest, CheckRemovePTMsWithModifications)
{
    // Remove these
    ASSERT_EQ("A.BCD.E", PSMDescription::removePTMs("A.n[1]BCD.E"));
    ASSERT_EQ("A.BCD.E", PSMDescription::removePTMs("A.BCDc[2].E"));
    // But not these
    ASSERT_EQ("A.nBCD.E", PSMDescription::removePTMs("A.nBCD.E"));
    ASSERT_EQ("A.BCcD.E", PSMDescription::removePTMs("A.BCc[2]D.E"));
    ASSERT_EQ("A.cBCDn.E", PSMDescription::removePTMs("A.c[1]BCDn[2].E"));
}
