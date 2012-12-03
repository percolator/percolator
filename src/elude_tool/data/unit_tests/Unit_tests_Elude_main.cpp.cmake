/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
/*
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/*
 * Main file for Google Test
 * Just run all tests
 */
#include "UnitTest_Elude_RetentionFeatures.cpp"
#include "UnitTest_Elude_EludeCaller.cpp" 
#include "UnitTest_Elude_DataManager.cpp"
#include "UnitTest_Elude_LibSVRModel.cpp"
#include "UnitTest_Elude_RetentionModel.cpp"


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
