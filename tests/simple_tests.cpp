#include <gtest/gtest.h>
#include "MassHandler.h"
#include "Normalizer.h"

namespace {

class MassHandlerTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      MassHandler::setMonoisotopicMass(true);
    }

    virtual void TearDown() {}
};

TEST_F(MassHandlerTest, BorderTest)
{
  EXPECT_EQ(0, MassHandler::massDiff(0.0, 0.0, 1, "ABCDEF"));
}

TEST_F(MassHandlerTest, AnotherBorderTest) {
  EXPECT_EQ(1, MassHandler::massDiff(0.0, 0.0, 1, "ABCDEF"));
}

class NormalizerTest : public ::testing::Test {
  protected:
  virtual void SetUp() {
    Normalizer::setType(Normalizer::UNI);
    cout << ";)" << endl;
    n = Normalizer::getNormalizer();
  }

  virtual void TearDown() {delete n;}
  Normalizer *n;
};

TEST_F(NormalizerTest, BorderTest) {
  vector<double> in, out, out2;
  for(int i = 0; i < 100; i++) {
    in.push_back(i / 2.0);
    cout << in[i] << " ";
  }

  n->normalizeweight(in, out);
  n->unnormalizeweight(out, out2);
  ASSERT_TRUE(in == out2);
}

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
