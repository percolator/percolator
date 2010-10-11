#include <gtest/gtest.h>
#include "MassHandler.h"
#include "Normalizer.h"

//using namespace std;

class MassHandlerTest : public ::testing::Test {
     protected:
          virtual void SetUp() {
               MassHandler::setMonoisotopicMass(true);
          }

          virtual void TearDown() {}
};

TEST_F(MassHandlerTest, BorderTest) {
     EXPECT_EQ(0, 0);
}

class NormalizerTest : public ::testing::Test {
     protected:
          virtual void SetUp() {
               Normalizer::setType(Normalizer::UNI);
               cout << ";)" << endl;
               n = Normalizer::getNormalizer();
          }

          virtual void TearDown() {
               delete n;
          }
          Normalizer* n;
};


int main(int argc, char** argv) {
     ::testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}  
