#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>    // for std::iota
#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse> 
#include <chrono>
#include "Globals.h"
#include "IsotonicPEP.h"


/**
 * Perform standard PAVA (pool-adjacent-violators) to enforce a non-decreasing sequence
 * on the input 'values'.  
 * Here we extended the method to also merge block if average is outside 
 * the range [min_value, max_value] 
 *
 * Each final "block" is repeated block.count times in the output.
 *
 * Typically used in 1D for y-values only.
 */
std::vector<double> PavaRegression::pavaNonDecreasingRanged(
    const std::vector<double>& values,
    const double min_value,
    const double max_value
) const
{
    const int n = static_cast<int>(values.size());
    if (n == 0) {
        return {};
    }

    // We'll store blocks on a stack
    struct MergeBlock {
        double sum;
        int    count;
        double avg;
    };

    std::vector<MergeBlock> stack;
    stack.reserve(n);

    // 1. Left to right
    for (int i = 0; i < n; ++i) {
        MergeBlock newBlock { values[i], 1, values[i] };
        stack.push_back(newBlock);

        // 2. Merge while there's a violation of non-decreasing
        while (stack.size() > 1) {
            auto &top    = stack.back();
            auto &secTop = stack[stack.size() - 2];
            double mergedSum   = secTop.sum + top.sum;
            int    mergedCount = secTop.count + top.count;
            double mergedAvg   = mergedSum / mergedCount;

            // if (( secTop.avg > top.avg) || ( mergedAvg < min_value ) || ( mergedAvg > max_value )) {
            if ( secTop.avg > top.avg) {
                stack.pop_back();
                stack.pop_back();
                stack.push_back({ mergedSum, mergedCount, mergedAvg });
            } else {
                break;
            }
        }
    }

    // 3. Expand final solution
    std::vector<double> result;
    result.reserve(n); // approximate; actual size = sum(counts)

    for (auto &b : stack) {
        double val = std::max( std::min(b.avg, max_value), min_value);
        for (int c = 0; c < b.count; ++c) {
            result.push_back(val);
        }
    }
    return result;
}

BinnedData IsplineRegression::bin_data_weighted(const std::vector<double>& x, const std::vector<double>& y, int max_bins) const {
    size_t n = x.size();
    if (n <= static_cast<size_t>(max_bins)) {
        return {x, y, std::vector<double>(n, 1.0)};
    }

    // Step 1: Smooth y using moving average to estimate local P(y=1)
    std::vector<double> p(n, 0.0);
    const size_t window = std::max<size_t>(10, n / 200); // around 0.5% of data
    for (size_t i = 0; i < n; ++i) {
        size_t start = (i >= window/2) ? i - window/2 : 0;
        size_t end = std::min(n, i + window/2);
        double sum = std::accumulate(y.begin() + start, y.begin() + end, 0.0);
        p[i] = sum / (end - start);
    }

    // Step 2: Compute cumulative *change* in estimated P(y=1)
    std::vector<double> change(n, 0.0);
    for (size_t i = 1; i < n; ++i) {
        change[i] = change[i - 1] + std::abs(p[i] - p[i - 1]);
    }
    double total_change = change.back();
    if (total_change == 0.0) total_change = 1.0;

    // Step 3: Create bins where change crosses thresholds
    std::vector<double> x_binned, y_binned, weights;
    size_t bin_start = 0;
    size_t target_bin = 1;

    for (size_t i = 1; i < n; ++i) {
        if (change[i] >= (target_bin * total_change / max_bins)) {
            size_t count = i - bin_start;
            if (count > 0) {
                double x_sum = std::accumulate(x.begin() + bin_start, x.begin() + i, 0.0);
                double y_sum = std::accumulate(y.begin() + bin_start, y.begin() + i, 0.0);
                x_binned.push_back(x_sum / count);
                y_binned.push_back(y_sum / count);
                weights.push_back(static_cast<double>(count));
            }
            bin_start = i;
            target_bin++;
        }
    }

    // Final bin
    if (bin_start < n) {
        size_t count = n - bin_start;
        double x_sum = std::accumulate(x.begin() + bin_start, x.end(), 0.0);
        double y_sum = std::accumulate(y.begin() + bin_start, y.end(), 0.0);
        x_binned.push_back(x_sum / count);
        y_binned.push_back(y_sum / count);
        weights.push_back(static_cast<double>(count));
    }

    return {x_binned, y_binned, weights};
}

BinnedData old_bin_data_weighted(
    const std::vector<double>& x, const std::vector<double>& y, int max_bins)
{
    size_t n = x.size();
    if (n <= static_cast<size_t>(max_bins)) {
        return {x, y, std::vector<double>(x.size(), 1.0)}; // No binning needed
    }

    std::vector<double> x_binned(max_bins, 0.0);
    std::vector<double> y_binned(max_bins, 0.0);
    std::vector<double> weights(max_bins, 0.0);

    double bin_size = static_cast<double>(n) / max_bins;

    for (size_t i = 0; i < n; ++i) {
        int bin = std::min(static_cast<int>(i / bin_size), max_bins - 1);
        x_binned[bin] += x[i];
        y_binned[bin] += y[i];
        weights[bin] += 1.0;
    }

    // Finalize average per bin
    std::vector<double> x_avg, y_avg, w_out;
    for (int i = 0; i < max_bins; ++i) {
        if (weights[i] > 0) {
            x_avg.push_back(x_binned[i] / weights[i]);
            y_avg.push_back(y_binned[i] / weights[i]);
            w_out.push_back(weights[i]);
        }
    }

    return {x_avg, y_avg, w_out};
}

std::vector<double> IsplineRegression::fit_y(const std::vector<double>& y,
    double min_val,
    double max_val) const {

    if (VERB > 2)
        std::cerr << "[TIMING] entering fit_y\n";

    int n = y.size();
    std::vector<double> x(n);
    std::iota(x.begin(), x.end(), 0);
    return fit_xy(x, y, min_val, max_val);
}

std::vector<double> IsplineRegression::fit_xy(const std::vector<double>& x, const std::vector<double>& y,
                                              double min_val, double max_val) const {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    if (VERB > 2) std::cerr << "[TIMING] entering sparse fit_xy\n";

    assert(x.size() == y.size());

    BinnedData tmp = bin_data_weighted(x, y);
    auto x_used = tmp.x;
    auto y_used = tmp.y;
    auto weights = tmp.weights;

    if (x_used.empty() || y_used.empty()) return {};

    int n = x_used.size();
    if (n == 0) return {};

    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    double xscale = xmax - xmin;

    std::vector<double> x_norm = normalize_vector(x_used, xmin, xscale);

    int num_knots = std::min(static_cast<int>(std::sqrt(n)), 50);
    std::vector<double> knots = compute_knots(x_norm, y_used, num_knots);

    if (VERB > 3) {
        std::cerr << "[DEBUG] Knots (" << knots.size() << "): ";
        for (double k : knots) std::cerr << k << " ";
        std::cerr << "\n";
    }


    Eigen::SparseMatrix<double, Eigen::RowMajor> X = build_design_matrix(x_norm, knots, n, num_knots);
    fill_intercept(X);

    if (VERB > 3) {
        std::cerr << "[DEBUG] First 50 rows of design matrix (non-zero entries):\n";
        for (int i = 0; i < std::min(50, n); ++i) {
            std::cerr << "Row " << i << ": ";
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(X, i); it; ++it)
                std::cerr << "(" << it.col() << ", " << it.value() << ") ";
            std::cerr << "\n";
        }
    }


    Eigen::VectorXd yvec = Eigen::Map<const Eigen::VectorXd>(y_used.data(), y_used.size());
    Eigen::VectorXd sqrt_weights_vec = Eigen::Map<const Eigen::VectorXd>(weights.data(), n).array().sqrt();

    apply_weights(X, yvec, sqrt_weights_vec);

    std::vector<bool> constrained(num_knots + 1, true);
    constrained[0] = false;

    auto box_start = high_resolution_clock::now();
    Eigen::VectorXd coeffs = box_lsq_ldlt(X, yvec, constrained);
    auto box_end = high_resolution_clock::now();
    if (VERB > 3) {
        std::cerr << "[DEBUG] Coefficients:\n";
        for (int i = 0; i < coeffs.size(); ++i)
            std::cerr << "  coeff[" << i << "] = " << coeffs[i] << "\n";
    }
    if (VERB > 2)
        std::cerr << "[TIMING] box_lsq_ldlt sparse duration: " << std::chrono::duration<double>(box_end - box_start).count() << " seconds\n";

    int n_full = x.size();
    std::vector<double> x_full_norm = normalize_vector(x, xmin, xscale);

    Eigen::MatrixXd X_full_dense = build_dense_matrix(x_full_norm, knots, n_full, num_knots);

    Eigen::VectorXd yfit_full = X_full_dense * coeffs;

    if (VERB > 3) {
        std::cerr << "[DEBUG] Binned y values (first 20): ";
        for (size_t i = 0; i < std::min(y_used.size(), size_t(20)); ++i)
            std::cerr << y_used[i] << " ";
        std::cerr << "\n";

        std::cerr << "[DEBUG] First 10 predictions (before clamp):\n";
        for (int i = 0; i < std::min(10, n_full); ++i)
            std::cerr << "  yfit[" << i << "] = " << yfit_full[i] << "\n";
    }

    std::vector<double> result(n_full);
    std::transform(yfit_full.data(), yfit_full.data() + n_full, result.begin(), [&](double v) {
        return clamp(v, min_val, max_val);
    });

    for (size_t i = 1; i < result.size(); ++i)
        result[i] = std::max(result[i], result[i - 1]);

    auto end = high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] sparse fit_xy total duration: " << std::chrono::duration<double>(end - start).count() << " seconds\n";

    return result;
}

std::vector<double> IsplineRegression::normalize_vector(const std::vector<double>& x, double xmin, double xscale) const {
    std::vector<double> norm(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        norm[i] = (x[i] - xmin) / xscale;
    return norm;
}

std::vector<double> IsplineRegression::compute_knots(const std::vector<double>& /* x_norm */,
                                                     const std::vector<double>& /* y */,
                                                     int num_knots) const {
    std::vector<double> knots(num_knots + 1);
    for (int i = 0; i <= num_knots; ++i)
    knots[i] = static_cast<double>(i) / num_knots;
    return knots;
}

Eigen::SparseMatrix<double, Eigen::RowMajor> IsplineRegression::build_design_matrix(const std::vector<double>& x_norm, const std::vector<double>& knots, int n, int num_knots) const {
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_per_thread(omp_get_max_threads());
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& local_triplets = triplets_per_thread[tid];
        local_triplets.reserve(n / omp_get_max_threads() + 8);

        #pragma omp for nowait
        for (int i = 0; i < n; ++i) {
            double xi = x_norm[i];
            auto it = std::upper_bound(knots.begin(), knots.end(), xi);
            int j = std::max(0, static_cast<int>(it - knots.begin()) - 1);
            if (j < num_knots) {
                double val = cubic_ispline(xi, knots[j], knots[j + 1]);
                if (val > 0.0)
                    local_triplets.emplace_back(i, j + 1, val);
            }
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    for (const auto& vec : triplets_per_thread)
        triplets.insert(triplets.end(), vec.begin(), vec.end());

    Eigen::SparseMatrix<double, Eigen::RowMajor> X(n, num_knots + 1);
    X.reserve(Eigen::VectorXi::Constant(n, 2));
    X.setFromTriplets(triplets.begin(), triplets.end());
    X.makeCompressed();
    return X;
}

void IsplineRegression::fill_intercept(Eigen::SparseMatrix<double, Eigen::RowMajor>& X) const {
    for (int i = 0; i < X.rows(); ++i)
        X.coeffRef(i, 0) = 1.0;
}

void IsplineRegression::apply_weights(Eigen::SparseMatrix<double, Eigen::RowMajor>& X, Eigen::VectorXd& yvec, const Eigen::VectorXd& sqrt_weights_vec) const {
    for (int i = 0; i < X.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(X, i); it; ++it) {
            it.valueRef() *= sqrt_weights_vec[it.row()];
        }
    }
    yvec = yvec.cwiseProduct(sqrt_weights_vec);
}

Eigen::MatrixXd IsplineRegression::build_dense_matrix(const std::vector<double>& x_full_norm, const std::vector<double>& knots, int n_full, int num_knots) const {
    Eigen::MatrixXd X_full_dense = Eigen::MatrixXd::Zero(n_full, num_knots + 1);
    for (int i = 0; i < n_full; ++i)
        X_full_dense(i, 0) = 1.0;

    for (int i = 0; i < n_full; ++i) {
        double xi = x_full_norm[i];
        int j = 0;
        while (j + 1 < static_cast<int>(knots.size()) && xi > knots[j + 1])
            ++j;

        int width = 5;
        for (int k = std::max(0, j - width); k <= std::min(j + width, num_knots - 1); ++k) {
            double val = cubic_ispline(xi, knots[k], knots[k + 1]);
            if (val > 0.0)
                X_full_dense(i, k + 1) = val;
        }
    }
    return X_full_dense;
}

double IsplineRegression::cubic_ispline(double x, double left, double right) const {
    if (x < left) return 0.0;
    if (x >= right) return 1.0;
    double u = (x - left) / (right - left);
    return 3 * u * u - 2 * u * u * u;
}

Eigen::VectorXd IsplineRegression::box_lsq_ldlt(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const std::vector<bool>& constrained,
    double lambda,
    int max_iter,
    double tol
) const {
    using namespace Eigen;

    const int n = X.rows();
    const int p = X.cols();

    VectorXd x = VectorXd::Zero(p);
    x[0] = 0.5; // Warm-start the intercept
    VectorXd w = X.transpose() * (y - X * x);

    std::vector<bool> passive(p, false);
    std::vector<int> passive_set;

    // Pre-activate intercept and first basis function
    passive[0] = true;
    passive[1] = true;
    passive_set.push_back(0);
    passive_set.push_back(1);

    for (int iter = 0; iter < max_iter; ++iter) {
        int idx = -1;
        double max_w = 0;
        for (int i = 0; i < p; ++i) {
            if (!passive[i] && constrained[i] && w[i] > max_w) {
                max_w = w[i];
                idx = i;
            }
        }

        if (idx == -1)
            break;

        passive[idx] = true;
        passive_set.push_back(idx);

        bool inner_loop = true;
        while (inner_loop) {
            MatrixXd X_p(n, passive_set.size());
            for (size_t i = 0; i < passive_set.size(); ++i)
                X_p.col(i) = X.col(passive_set[i]);

            VectorXd z_p = X_p.transpose() * y;
            MatrixXd XtX = X_p.transpose() * X_p;
            if (lambda > 0) XtX += lambda * MatrixXd::Identity(XtX.rows(), XtX.cols());

            VectorXd b_p = XtX.ldlt().solve(z_p);

            VectorXd z = VectorXd::Zero(p);
            for (size_t i = 0; i < passive_set.size(); ++i)
                z[passive_set[i]] = b_p[i];

            bool all_non_negative = true;
            for (int j : passive_set) {
                if (constrained[j] && z[j] < 0) {
                    all_non_negative = false;
                    break;
                }
            }

            if (all_non_negative) {
                x = z;
                inner_loop = false;
            } else {
                double alpha = 1.0;
                for (int j : passive_set) {
                    if (constrained[j] && z[j] < 0)
                        alpha = std::min(alpha, x[j] / (x[j] - z[j]));
                }
                x = x + alpha * (z - x);

                std::vector<int> new_passive;
                for (int j : passive_set) {
                    if (x[j] > tol) {
                        new_passive.push_back(j);
                    } else {
                        passive[j] = false;
                    }
                }
                passive_set = std::move(new_passive);
            }
        }

        w = X.transpose() * (y - X * x);
    }

    return x;
}

namespace util {
    template <typename T>
    const T& clamp(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

AdaptiveIsplineRegression::AdaptiveIsplineRegression(int degree) : degree_(degree) {}

std::vector<double> AdaptiveIsplineRegression::fit_y(const std::vector<double>& y, double min_val, double max_val) const {
    std::vector<double> x(y.size());
    std::iota(x.begin(), x.end(), 0.0);
    return fit_xy(x, y, min_val, max_val);
}

AdaptiveIsplineRegression::BinnedData AdaptiveIsplineRegression::adaptive_bin(const std::vector<double>& x, const std::vector<double>& y, int max_bins) const {
    size_t n = x.size();
    std::vector<double> x_binned, y_binned, weights;

    size_t bin_start = 0;
    for (size_t i = 1; i <= n; ++i) {
        if (i == n || (y[i] != y[bin_start] && (i - bin_start >= n / max_bins || y[i] == 0))) {
            double x_avg = std::accumulate(x.begin() + bin_start, x.begin() + i, 0.0) / (i - bin_start);
            double y_avg = std::accumulate(y.begin() + bin_start, y.begin() + i, 0.0) / (i - bin_start);
            if (VERB > 2) std::cerr << "Creating bin from " << bin_start << " to " << i - 1 << ", x_avg: " << x_avg << ", y_avg: " << y_avg << "\n";
            x_binned.push_back(x_avg);
            y_binned.push_back(y_avg);
            weights.push_back(i - bin_start);
            bin_start = i;
        }
    }
    return {x_binned, y_binned, weights};
}

std::vector<double> AdaptiveIsplineRegression::compute_adaptive_knots(const std::vector<double>& x, const std::vector<double>& /* y */, int num_knots) const {
    std::vector<double> knots;
    knots.push_back(x.front());
    for (int i = 1; i < num_knots; ++i) {
        double q = static_cast<double>(i) / num_knots;
        size_t idx = q * (x.size() - 1);
        knots.push_back(x[idx]);
    }
    knots.push_back(x.back());
    if (VERB > 1) {
        std::cerr << "Knots: ";
        for (const auto& k : knots) std::cerr << k << " ";
        std::cerr << "\n";
    }
    return knots;
}

Eigen::VectorXd AdaptiveIsplineRegression::fit_spline(const BinnedData& data, const std::vector<double>& knots, double lambda) const {
    int n = data.x.size(), k = knots.size();
    Eigen::MatrixXd X(n, k);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < k - 1; ++j)
            X(i, j) = quadratic_ispline(data.x[i], knots[j], knots[j + 1]);
    X.col(k - 1).setOnes();

    Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(data.y.data(), n);
    Eigen::VectorXd w = Eigen::Map<const Eigen::VectorXd>(data.weights.data(), n);
    Eigen::VectorXd coeffs = (X.transpose() * w.asDiagonal() * X + lambda * Eigen::MatrixXd::Identity(k, k)).ldlt().solve(X.transpose() * w.asDiagonal() * y);

    if (VERB > 1) std::cerr << "Spline coefficients: " << coeffs.transpose() << "\n";

    return coeffs;
}

double AdaptiveIsplineRegression::quadratic_ispline(double x, double left, double right) const {
    if (x < left) return 0.0;
    if (x > right) return 1.0;
    double u = (x - left) / (right - left);
    return u * u;
}

std::vector<double> AdaptiveIsplineRegression::fit_xy(const std::vector<double>& x, const std::vector<double>& y, double min_val, double max_val) const {
    assert(x.size() == y.size());

    auto data = adaptive_bin(x, y, 2500);
    auto knots = compute_adaptive_knots(data.x, data.y, std::min(50, (int)std::sqrt(data.x.size())));

    Eigen::VectorXd coeffs = fit_spline(data, knots, 1e-8);

    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        double pred = coeffs.tail(1)(0);
        for (size_t j = 0; j < knots.size() - 1; ++j)
            pred += coeffs(j) * quadratic_ispline(x[i], knots[j], knots[j + 1]);
        result[i] = util::clamp(pred, min_val, max_val);
    }

    return PavaRegression().fit_xy(x, result, min_val, max_val);
}

InferPEP::InferPEP(bool use_ispline)
{
    if (use_ispline) {
        regressor_ptr_ = std::make_unique<AdaptiveIsplineRegression>(2);
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using I-Splines" << std::endl;
        }                
    } else {
        regressor_ptr_ = std::make_unique<PavaRegression>();
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using PAVA" << std::endl;
        }
    }
}

std::vector<double> InferPEP::q_to_pep(const std::vector<double>& q_values) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }
    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }   
    return regressor_ptr_->fit_y(raw_pep);
}

std::vector<double> InferPEP::qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }

    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }

    return regressor_ptr_->fit_xy(scores, raw_pep);
}

std::vector<double> InferPEP::tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores) {

    using namespace std::chrono;
    auto start = std::chrono::high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] entering tdc_to_pep\n";

    double epsilon = 1e-20;
    auto is_dec = is_decoy;
    is_dec.insert(is_dec.begin(), 0.5);

    std::vector<double> decoy_rate;
    if (!scores.empty()) {
        if (VERB > 2)
            std::cerr << "[TIMING] choosing fit_xy\n";
        auto sc = scores;
        sc.insert(sc.begin(), sc[0]);
        decoy_rate = regressor_ptr_->fit_xy(sc, is_dec);
    } else {
        if (VERB > 2)
            std::cerr << "[TIMING] choosing fit_y\n";
        decoy_rate = regressor_ptr_->fit_y(is_dec);
    }
    decoy_rate.erase(decoy_rate.begin());

    std::vector<double> pep_iso;
    // auto i_d = is_decoy.begin(); double tar = 0.; double dec = 0.; double errors = 0.;
    for (auto& dp : decoy_rate) {
        if (dp > 1. - epsilon)
            dp = 1. - epsilon;
        double pep = dp / (1 - dp);
        if (pep > 1.)
            pep = 1.;
        pep_iso.push_back(pep);
        // if (*i_d>0.5) {dec+=1.;} else {tar+=1.; errors += pep;} i_d++;
        // cerr << dp << "\t" << pep << "\t" << tar << "\t" << dec << "\t" << errors << "\t" << (dec+0.5)/(tar) << "\t" << errors/(tar) << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    if (VERB > 2)
        std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";

    return pep_iso;
}
