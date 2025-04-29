#include <gtest/gtest.h>
#include "IsotonicPEP.h" // Your header file with IsplineRegression

class IsplineRegressionTest : public ::testing::Test {
protected:
    IsplineRegression model;

    void SetUp() override {
        model = IsplineRegression(); // default degree=2
    }
};

TEST_F(IsplineRegressionTest, MonotonicityFitY) {
    // Test that fit_y returns a non-decreasing sequence
    std::vector<double> y = {0.1, 0.05, 0.2, 0.15, 0.3, 0.25, 0.5};
    std::vector<double> fitted = model.fit_y(y, 0.0, 1.0);

    ASSERT_EQ(fitted.size(), y.size());

    for (size_t i = 1; i < fitted.size(); ++i) {
        ASSERT_GE(fitted[i], fitted[i-1]) << "Violation of monotonicity at position " << i;
    }
}

TEST_F(IsplineRegressionTest, FitXYBinaryInputs) {
    // Test behavior on a typical binary dataset (0s and 1s)
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 6};
    std::vector<double> y = {0, 0, 0, 1, 1, 1, 1};
    std::vector<double> fitted = model.fit_xy(x, y, 0.0, 1.0);

    ASSERT_EQ(fitted.size(), x.size());

    for (size_t i = 1; i < fitted.size(); ++i) {
        ASSERT_GE(fitted[i], fitted[i-1]) << "Violation of monotonicity at position " << i;
    }

    EXPECT_GE(fitted.front(), 0.0);
    EXPECT_LE(fitted.back(), 1.0);
}

TEST_F(IsplineRegressionTest, FitXYBinaryInputsLarge) {
    // Create a larger binary dataset: 1000 points
    const int N = 1000;
    std::vector<double> x(N);
    std::vector<double> y(N);

    for (int i = 0; i < N; ++i) {
        x[i] = i;
        if (i < N / 3)
            y[i] = 0; // early part 0s
        else if (i > 2 * N / 3)
            y[i] = 1; // late part 1s
        else
            y[i] = (i % 2 == 0) ? 0 : 1; // middle part, flip between 0 and 1 (forcing model to learn slope)
    }

    std::vector<double> fitted = model.fit_xy(x, y, 0.0, 1.0);

    ASSERT_EQ(fitted.size(), x.size());

    for (size_t i = 1; i < fitted.size(); ++i) {
        ASSERT_GE(fitted[i], fitted[i-1]) << "Violation of monotonicity at position " << i;
    }

    for (size_t i = 0; i < fitted.size(); ++i) {
        EXPECT_GE(fitted[i], 0.0) << "Negative fitted value at position " << i;
        EXPECT_LE(fitted[i], 1.0) << "Fitted value >1 at position " << i;
    }
}

TEST_F(IsplineRegressionTest, HandlesConstantInput) {
    // Test when all y are identical
    std::vector<double> y(100, 0.5);
    std::vector<double> fitted = model.fit_y(y);

    ASSERT_EQ(fitted.size(), y.size());

    for (size_t i = 0; i < fitted.size(); ++i) {
        ASSERT_NEAR(fitted[i], 0.5, 1e-6);
    }
}

TEST_F(IsplineRegressionTest, SmallDatasetUniformKnots) {
    // Dataset smaller than 200 points should trigger uniform spacing
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {0.1, 0.2, 0.2, 0.3, 0.4};
    std::vector<double> fitted = model.fit_xy(x, y);

    ASSERT_EQ(fitted.size(), x.size());
    for (size_t i = 1; i < fitted.size(); ++i) {
        ASSERT_GE(fitted[i], fitted[i-1]);
    }
}

TEST_F(IsplineRegressionTest, LargeDatasetAdaptiveKnots) {
    // Large dataset should trigger adaptive spacing
    std::vector<double> x(500), y(500);
    for (int i = 0; i < 500; ++i) {
        x[i] = i;
        y[i] = (i < 300) ? 0 : 1;
    }

    std::vector<double> fitted = model.fit_xy(x, y);

    ASSERT_EQ(fitted.size(), x.size());
    for (size_t i = 1; i < fitted.size(); ++i) {
        ASSERT_GE(fitted[i], fitted[i-1]);
    }
}