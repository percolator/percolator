#include <gtest/gtest.h>
#include <algorithm>
#include <iostream> // for std::cerr, std::endl
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include <tuple>
#include <string>

#include "IsotonicPEP.h"          // your class (InferPEP)
#include "MonotoneRegressor.h"    // RegressorType, MonotoneParams factory

namespace {

// Generate synthetic data as specified:
// - Top 1000 highest scores: all targets (non-decoy)
// - Next 1000: decoy prob ramps linearly 0 -> 0.5
// - Last 1000: decoy prob = 0.5
struct SynthData {
  std::vector<double> scores;   // larger is better (descending order)
  std::vector<double> is_decoy; // 1 for decoy, 0 for target
};

SynthData make_synth_data(unsigned seed = 1337) {
  const int n_block = 1000;
  const int N = 3 * n_block;

  // Scores: strictly decreasing, top are best
  std::vector<double> scores(N);
  for (int i = 0; i < N; ++i) scores[i] = double(N - i);

  std::mt19937 rng(seed);
  std::vector<double> is_decoy; is_decoy.reserve(N);

  // Block 1: prob=0 (all targets)
  for (int i = 0; i < n_block; ++i) is_decoy.push_back(0.0);

  // Block 2: prob ramps 0 -> 0.5
  for (int i = 0; i < n_block; ++i) {
    double p = 0.5 * (double(i) / (n_block - 1));
    std::bernoulli_distribution bern(p);
    is_decoy.push_back(bern(rng) ? 1.0 : 0.0);
  }

  // Block 3: prob = 0.5
  std::bernoulli_distribution bern_half(0.5);
  for (int i = 0; i < n_block; ++i)
    is_decoy.push_back(bern_half(rng) ? 1.0 : 0.0);

  return {scores, is_decoy};
}

SynthData make_sparse_top_decoy_data(unsigned seed = 2025) {
  const int n_top = 400;
  const int n_ramp = 1200;
  const int n_tail = 1200;
  const int N = n_top + n_ramp + n_tail;

  std::vector<double> scores(N);
  for (int i = 0; i < N; ++i) scores[i] = double(N - i);

  std::mt19937 rng(seed);
  std::vector<double> is_decoy; is_decoy.reserve(N);

  for (int i = 0; i < n_top; ++i) is_decoy.push_back(0.0);

  for (int i = 0; i < n_ramp; ++i) {
    double p = 0.20 * (double(i) / std::max(1, n_ramp - 1));
    std::bernoulli_distribution bern(p);
    is_decoy.push_back(bern(rng) ? 1.0 : 0.0);
  }

  std::bernoulli_distribution bern_tail(0.5);
  for (int i = 0; i < n_tail; ++i) {
    is_decoy.push_back(bern_tail(rng) ? 1.0 : 0.0);
  }

  return {scores, is_decoy};
}

// Given pep[i] for i=0..M-1 (aligned with scores/is_decoy),
// evaluate calibration at each prefix i where the i-th item is a TARGET:
//   fdr_hat(i)  = (# decoys up to i) / (# targets up to i)
//   pep_mean(i) = (sum of pep_j for targets j <= i) / (# targets up to i)
// We assert relative deviation <= rel_tol after a warmup (to reduce noise).
void check_calibration_q(const std::vector<double>& pep,
                         const std::vector<double>& is_decoy,
                         double rel_tol = 0.20)
{
  const size_t N = pep.size();

  // prefix counts
  std::vector<int> D(N,0), T(N,0);
  for (size_t i=0;i<N;++i) {
    D[i] = (i?D[i-1]:0) + (is_decoy[i] > 0.5 ? 1:0);
    T[i] = (i?T[i-1]:0) + (is_decoy[i] > 0.5 ? 0:1);
  }

  // raw fdrs with +0.5 offset
  std::vector<double> fdr(N,0.0);
  for (size_t i=0;i<N;++i) {
    fdr[i] = (T[i]>0) ? ( (double)D[i] + 0.5 ) / (double)T[i] : 0.0;
    if (fdr[i] > 1.0) fdr[i] = 1.0;
  }

  // backward smoothing into q-values
  std::vector<double> q(N,1.0);
  double run_min = 1.0;
  for (ptrdiff_t i=(ptrdiff_t)N-1; i>=0; --i) {
    run_min = std::min(run_min, fdr[(size_t)i]);
    q[(size_t)i] = run_min;
  }

  // compare mean PEP vs q
  constexpr int kWarmupTargets = 25;
  constexpr double kAbsTol = 2e-3;

  int targets=0; double pep_sum=0.0;
  for (size_t i=0;i<N;++i) {
    if (is_decoy[i] > 0.5) continue;
    pep_sum += pep[i];
    targets++;

    if (targets < kWarmupTargets) {
      continue;
    }

    double pep_mean = pep_sum / (double)targets;
    double qi = q[i];
    double denom = std::max(1e-6, std::min(pep_mean,qi));
    double rel_err = std::abs(pep_mean - qi)/denom;
    double abs_err = std::abs(pep_mean - qi);
    ASSERT_TRUE(rel_err < rel_tol || abs_err < kAbsTol)
        << "Poor calibration at index " << i
        << " pep_mean=" << pep_mean
        << " vs q[i]=" << qi
        << " rel_err=" << rel_err
        << " abs_err=" << abs_err;
  }

  EXPECT_NEAR(q.back(),
              pep_sum/std::max(1,targets),
              0.15); // soft tail check
}
} // namespace

// ------------------------------------------------------------------
// Parameterization: (use_ispline, use_fit_xy)
//   false = PAVA, true = ISPLINE
//   false = fit_y, true = fit_xy
// ------------------------------------------------------------------
class TdcToPepCalibrationTest
    : public ::testing::TestWithParam<std::tuple<bool, bool>> {};

// Custom parameter name generator for friendly test names
static std::string ParamNameGen(
    const ::testing::TestParamInfo<std::tuple<bool, bool>>& info) {
  const auto [use_ispline, use_fit_xy] = info.param;
  const char* reg = use_ispline ? "ISpline" : "Pava";
  const char* fit = use_fit_xy   ? "FitXY"   : "FitY";
  return std::string(reg) + "_" + fit;
}

INSTANTIATE_TEST_SUITE_P(
    AllMonoRegressors,
    TdcToPepCalibrationTest,
    ::testing::Combine(::testing::Bool(), ::testing::Bool()),
    ParamNameGen);

// --- Tests ---------------------------------------------------------

TEST_P(TdcToPepCalibrationTest, SyntheticRamp_CalibratedPEPs) {
  // 1) Synthesize
  auto data = make_synth_data(/*seed=*/1337);
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  // 2) Fit using chosen method

  std::vector<double> pep;
  if (use_fit_xy) {
    pep = infer.tdc_to_pep(data.is_decoy, data.scores);
  } else {
    pep = infer.tdc_to_pep(data.is_decoy);
  }

  // 4) Compare cumulative mean PEP among targets vs TDC FDR
  check_calibration_q(pep, data.is_decoy, /*rel_tol=*/3.0); // tighten as needed
}

TEST_P(TdcToPepCalibrationTest, GetParam_Sane_AfterFit) {
  auto data = make_synth_data(/*seed=*/42);
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  if (use_fit_xy) {
    infer.tdc_to_pep(data.is_decoy, data.scores);
  } else {
    infer.tdc_to_pep(data.is_decoy);
  }

  // Smoke: model remains usable
  auto pep = infer.tdc_to_pep(data.is_decoy, data.scores);
}

TEST_P(TdcToPepCalibrationTest, SparseHighScoreDecoys_DoNotForceLargeTailFloor) {
  auto data = make_sparse_top_decoy_data();
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  std::vector<double> pep = use_fit_xy
      ? infer.tdc_to_pep(data.is_decoy, data.scores)
      : infer.tdc_to_pep(data.is_decoy);

  double top_target_sum = 0.0;
  int top_target_count = 0;
  for (size_t i = 0; i < pep.size() && top_target_count < 100; ++i) {
    if (data.is_decoy[i] > 0.5) continue;
    top_target_sum += pep[i];
    ++top_target_count;
  }

  ASSERT_EQ(top_target_count, 100);
  EXPECT_LT(top_target_sum / top_target_count, 0.05);
}

// PAVA with tied x: observations with the same x are one block, so they
// get the same fitted value whatever their order in the input. They used
// to be fed one by one, in input order, and a run of equal x was split
// where its y values happened to be ordered.
TEST(PavaTiesTest, EqualXGetsOneValue) {
  auto pava = make_monotone_regressor(RegressorType::PAVA);
  // decreasing in x (the default): three observations at x = 3 with mean
  // 1/3, two at x = 2 with mean 1/2, one at x = 1
  const std::vector<double> x = {3.0, 3.0, 3.0, 2.0, 2.0, 1.0};
  const std::vector<double> y = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  const std::vector<double> fit = pava->fit_xy(x, y);
  ASSERT_EQ(6u, fit.size());
  // fed one by one this gave 0, 0.5, 0.5, 0.5, 0.5, 1
  const std::vector<double> expected = {1.0 / 3, 1.0 / 3, 1.0 / 3,
                                        0.5, 0.5, 1.0};
  for (size_t i = 0; i < fit.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected[i], fit[i]) << "observation " << i;
  }

  // a tie block that violates the order is merged with its neighbour as a
  // whole, weighted by its size: (1 + 1 + 0) at x = 3 and (0) at x = 2
  const std::vector<double> fit2 =
      pava->fit_xy({3.0, 3.0, 3.0, 2.0, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0});
  for (size_t i = 0; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(0.5, fit2[i]) << "observation " << i;
  }
  EXPECT_DOUBLE_EQ(1.0, fit2[4]);
}

TEST(PavaTiesTest, FitDoesNotDependOnTheOrderOfTiedObservations) {
  std::mt19937 rng(7);
  const int n = 2000;
  // scores on a coarse grid, so that most of them are tied; decoys are
  // more likely at low scores
  std::vector<double> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    const int level = static_cast<int>(rng() % 40);
    x[i] = 0.25 * level;
    y[i] = (rng() % 80) < static_cast<unsigned>(40 - level) ? 1.0 : 0.0;
  }
  auto pava = make_monotone_regressor(RegressorType::PAVA);
  const std::vector<double> fit = pava->fit_xy(x, y);
  // equal x, equal fit; and the fit does not increase with x
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int a, int b) { return x[a] > x[b]; });
  int distinct = 1;
  for (int k = 1; k < n; ++k) {
    if (x[order[k]] == x[order[k - 1]]) {
      EXPECT_EQ(fit[order[k - 1]], fit[order[k]]) << "x " << x[order[k]];
    } else {
      EXPECT_LE(fit[order[k - 1]], fit[order[k]]) << "x " << x[order[k]];
      ++distinct;
    }
  }
  EXPECT_EQ(40, distinct);
  // the fit keeps the sum of y
  EXPECT_NEAR(std::accumulate(y.begin(), y.end(), 0.0),
              std::accumulate(fit.begin(), fit.end(), 0.0), 1e-9);
  // another input order of the same observations gives the same fit
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  std::shuffle(perm.begin(), perm.end(), rng);
  std::vector<double> xp(n), yp(n);
  for (int i = 0; i < n; ++i) {
    xp[i] = x[perm[i]];
    yp[i] = y[perm[i]];
  }
  const std::vector<double> fitp = pava->fit_xy(xp, yp);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(fit[perm[i]], fitp[i]) << "observation " << perm[i];
  }
}

// Without ties the fit is what it was: that of the y values in the order
// of x, one block per observation.
TEST(PavaTiesTest, DistinctXIsTheFitOfTheSortedValues) {
  std::mt19937 rng(11);
  const int n = 500;
  std::vector<double> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    x[i] = n - i;  // descending, distinct
    y[i] = (rng() % n) < static_cast<unsigned>(i) ? 1.0 : 0.0;
  }
  auto pava = make_monotone_regressor(RegressorType::PAVA);
  EXPECT_EQ(pava->fit_y(y), pava->fit_xy(x, y));
}

// --tdc-pep --pava-pep: PSMs with the same score get the same PEP.
TEST(PavaTiesTest, TdcPepsOfEqualScoresAreEqual) {
  std::mt19937 rng(3);
  const int n = 3000;
  std::vector<double> scores(n), is_decoy(n);
  for (int i = 0; i < n; ++i) {
    const int level = 60 - i / 50;  // 50 PSMs per score, best first
    scores[i] = 0.5 * level;
    is_decoy[i] = (rng() % 120) < static_cast<unsigned>(60 - level) ? 1.0
                                                                     : 0.0;
  }
  InferPEP infer(/*use_ispline=*/false);
  const std::vector<double> pep = infer.tdc_to_pep(is_decoy, scores);
  ASSERT_EQ(static_cast<size_t>(n), pep.size());
  for (int i = 1; i < n; ++i) {
    if (scores[i] == scores[i - 1]) {
      EXPECT_EQ(pep[i - 1], pep[i]) << "PSM " << i;
    } else {
      EXPECT_LE(pep[i - 1], pep[i]) << "PSM " << i;
    }
  }
}
