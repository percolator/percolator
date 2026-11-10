#include <gtest/gtest.h>
#include "ISplineTRRRegressor.h"
#include "IsotonicPEP.h" // Your header file with IsplineRegression
#include "MyException.h"
#include <limits>
#include <string>

class IsplineRegressionTest : public ::testing::Test {
protected:
    ISplineTRRRegressor model;

    void SetUp() override {
        model = ISplineTRRRegressor(); // default degree=2
    }
};

TEST_F(IsplineRegressionTest, IsplineDesignMatrixMonotonicity) {
    // Verify that I-spline columns are monotone non-decreasing in x
    std::vector<double> scores;
    for (int i = 0; i <= 100; ++i) scores.push_back(i / 100.0);

    // Build a clamped knot vector: degree+1 boundary copies + internal knots
    std::vector<double> knots = {0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1};
    DenseMatrix X = build_ispline_design(scores, 3, knots, true);

    // Column 0 is intercept (all 1s); remaining columns are I-splines
    for (int col = 1; col < X.cols(); ++col) {
        for (int row = 1; row < X.rows(); ++row) {
            EXPECT_GE(X(row, col), X(row - 1, col) - 1e-12)
                << "I-spline column " << col << " not monotone at row " << row;
        }
        // Should start near 0 and end near 1
        EXPECT_NEAR(X(0, col), 0.0, 1e-10);
        EXPECT_NEAR(X(X.rows() - 1, col), 1.0, 1e-10);
    }
}

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
    std::vector<double> y(10000, 0.5);
    std::vector<double> fitted = model.fit_y(y);

    ASSERT_EQ(fitted.size(), y.size());
    /*
    for (size_t i = 0; i < fitted.size(); ++i) {
        std::cerr << i << "\t" << fitted[i] << std::endl;
    }
    */
    for (size_t i = 0; i < fitted.size(); ++i) {
        ASSERT_NEAR(fitted[i], 0.5, 3e-2);
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
// Scores of 1e308 and -1e308 have an infinite span: every fitted value was
// NaN, and qvality printed "-nan" PEPs with exit status 0.
TEST_F(IsplineRegressionTest, RejectsScoresWithAnInfiniteSpan) {
    const std::vector<double> y = {0.0, 0.0, 0.5, 1.0, 1.0};
    EXPECT_THROW(model.fit_xy({1e308, 3.0, 2.0, 1.0, -1e308}, y), MyException);
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(model.fit_xy({inf, 3.0, 2.0, 1.0, 0.0}, y), MyException);
    try {
        model.fit_xy({1e308, 3.0, 2.0, 1.0, -1e308}, y);
        FAIL() << "no exception";
    } catch (const MyException& e) {
        EXPECT_NE(std::string::npos,
                  std::string(e.what()).find("range from -1e+308 to 1e+308"))
            << e.what();
    }
    // the largest finite span is accepted
    std::vector<double> fitted;
    EXPECT_NO_THROW(fitted = model.fit_xy({8e307, 3.0, 2.0, 1.0, -8e307}, y));
    for (double v : fitted) EXPECT_TRUE(v >= 0.0 && v <= 1.0) << v;
}

// The solver's converged flag was not read: a fit that stopped at the
// iteration limit was used without a message.
TEST_F(IsplineRegressionTest, WarnsWhenTheSolverDoesNotConverge) {
    std::vector<double> x, y;
    for (int i = 0; i < 200; ++i) {
        x.push_back(200.0 - i);
        y.push_back(i < 100 ? 0.0 : 1.0);
    }
    testing::internal::CaptureStderr();
    const std::vector<double> converged = model.fit_xy(x, y);
    EXPECT_EQ("", testing::internal::GetCapturedStderr());

    // this fit needs more than two iterations
    MonotoneParams params = model.params();
    params.max_solver_iters = 2;
    model.set_params(params);
    testing::internal::CaptureStderr();
    const std::vector<double> stopped = model.fit_xy(x, y);
    const std::string warning = testing::internal::GetCapturedStderr();
    EXPECT_EQ("Warning: the I-spline regression of the PEPs did not converge "
              "(2 iterations, 18 coefficients); the PEPs may be inaccurate.\n",
              warning);
    EXPECT_NE(converged, stopped);
    // the stopped fit is still monotone and within the bounds
    for (size_t i = 0; i < stopped.size(); ++i) {
        EXPECT_TRUE(stopped[i] >= 0.0 && stopped[i] <= 1.0) << i;
        if (i > 0) EXPECT_GE(stopped[i], stopped[i - 1] - 1e-12) << i;
    }
}

// the solver reports whether it reached the optimum
TEST_F(IsplineRegressionTest, SolverReportsConvergence) {
    QuadraticModel q;
    q.G = DenseMatrix(2, 2);
    q.G(0, 0) = 2.0; q.G(0, 1) = 1.0; q.G(1, 0) = 1.0; q.G(1, 1) = 2.0;
    q.h = {1.0, -3.0};
    q.c = 0.0;
    const Vec lb(2, 0.0), ub(2, 10.0);
    const BoxQPResult solved = solveBoxQP(q, lb, ub, Vec(2, 5.0));
    EXPECT_TRUE(solved.converged);
    EXPECT_DOUBLE_EQ(0.5, solved.x[0]);
    EXPECT_DOUBLE_EQ(0.0, solved.x[1]);
    const BoxQPResult stopped = solveBoxQP(q, lb, ub, Vec(2, 5.0), 1);
    EXPECT_FALSE(stopped.converged);
    EXPECT_EQ(1, stopped.iters);
    // not positive definite
    q.G(1, 1) = -2.0;
    EXPECT_FALSE(solveBoxQP(q, lb, ub, Vec(2, 5.0)).converged);
}
