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
        return {x, y, std::vector<double>(x.size(), 1.0)}; // No binning needed, uniform weights
    }

    size_t bin_size = n / max_bins;
    std::vector<double> x_binned, y_binned, weights;
    x_binned.reserve(max_bins);
    y_binned.reserve(max_bins);
    weights.reserve(max_bins);

    for (size_t i = 0; i < n; i += bin_size) {
        size_t end = std::min(i + bin_size, n);
        double x_sum = 0.0, y_sum = 0.0;
        size_t count = end - i;
        for (size_t j = i; j < end; ++j) {
            x_sum += x[j];
            y_sum += y[j];
        }
        if (count > 0) {
            x_binned.push_back(x_sum / count);
            y_binned.push_back(y_sum / count);
            weights.push_back(static_cast<double>(count));
        }
    }

    return {x_binned, y_binned, weights};
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
    double min_val,
    double max_val) const {

    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] entering sparse fit_xy\n";

    assert(x.size() == y.size());

    BinnedData tmp = bin_data_weighted(x, y);
    auto x_used = tmp.x;
    auto y_used = tmp.y;
    auto weights = tmp.weights;

    if (x_used.empty() || y_used.empty()) {
        return {}; // empty input
    }

    int n = x_used.size();
    if (n == 0) return {};

    // Normalize x to [0,1]
    std::vector<double> x_norm(n);
    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    double xscale = xmax - xmin;
    for (int i = 0; i < n; ++i)
        x_norm[i] = (x[i] - xmin) / xscale;

    // Determine number of knots
    const int MAX_KNOTS = 50;
    int num_knots = std::min(static_cast<int>(std::sqrt(n)), MAX_KNOTS);

    std::vector<double> knots(num_knots + 1);

    if (n < 200) {
        // --- Small dataset: uniform knots ---
        for (int i = 0; i <= num_knots; ++i)
            knots[i] = static_cast<double>(i) / num_knots;
    } else {
        // --- Large dataset: adaptive knots based on local y movement ---

        // 1. Sort x_norm and y together
        std::vector<std::pair<double, double>> xy(n);
        for (int i = 0; i < n; ++i)
            xy[i] = {x_norm[i], y[i]};
        std::sort(xy.begin(), xy.end());

        // 2. Compute moving average of y
        std::vector<double> y_moving_avg(n, 0.0);
        const int window = std::max(5, n / 200); // window ~ 0.5% of data
        for (int i = 0; i < n; ++i) {
            int start = std::max(0, i - window/2);
            int end = std::min(n, i + window/2);
            double sum = 0;
            for (int j = start; j < end; ++j)
                sum += xy[j].second;
            y_moving_avg[i] = sum / (end - start);
        }

        // 3. Compute cumulative "action" score
        std::vector<double> action(n, 0.0);
        action[0] = 0.0;
        for (int i = 1; i < n; ++i)
            action[i] = action[i-1] + std::abs(y_moving_avg[i] - y_moving_avg[i-1]);

        // Normalize action to [0,1]
        double total_action = action.back();
        if (total_action == 0.0) total_action = 1.0; // avoid div-by-zero
        for (int i = 0; i < n; ++i)
            action[i] /= total_action;

        // 4. Distribute knots, with a couple ones fixed 
        knots[0] = 0.0;
        knots[1] = 0.02;  // very close to start
        knots[num_knots-1] = 0.98; // very close to end
        knots[num_knots] = 1.0;

        // 5. Fill adaptive inside
        int adaptive_knots = num_knots - 3; // we reserved 4 fixed, but there are n-1 intervalls
        for (int k = 1; k <= adaptive_knots; ++k) {
            double target = static_cast<double>(k) / (adaptive_knots + 1);
            auto it = std::lower_bound(action.begin(), action.end(), target);
            int idx = std::min(static_cast<int>(it - action.begin()), n-1);
            knots[k+1] = xy[idx].first; // Fill inside (shifted by +1)
        }
    }

    // Build Sparse Matrix of Cubic I-splines
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_per_thread(omp_get_max_threads());

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& local_triplets = triplets_per_thread[tid];
        local_triplets.reserve(n / omp_get_max_threads() + 8);
    
        #pragma omp for nowait
        for (int i = 0; i < n; ++i) {
            double xi = x_norm[i];
    
            // Locate interval
            auto it = std::upper_bound(knots.begin(), knots.end(), xi);
            int j = std::max(0, static_cast<int>(it - knots.begin()) - 1);        
            if (j < num_knots) {
                double val = cubic_ispline(xi, knots[j], knots[j+1]);
                if (val > 0.0)
                    local_triplets.emplace_back(i, j + 1, val);
            }
        }
    }
    
    // Merge all local vectors
    std::vector<Eigen::Triplet<double>> triplets;
    size_t total_size = 0;
    for (const auto& vec : triplets_per_thread)
        total_size += vec.size();
    
    triplets.reserve(total_size);
    
    for (const auto& vec : triplets_per_thread)
        triplets.insert(triplets.end(), vec.begin(), vec.end());    
    
    Eigen::SparseMatrix<double> X(n, num_knots + 1);
    X.reserve(Eigen::VectorXi::Constant(n, 2)); // Approx two non-zeros per row
    X.setFromTriplets(triplets.begin(), triplets.end());

    // Fill intercept column (first column, all ones)
    for (int i = 0; i < n; ++i)
        X.coeffRef(i, 0) = 1.0;

    Eigen::VectorXd yvec = Eigen::Map<const Eigen::VectorXd>(y_used.data(), y_used.size());

    // Reweight rows
    Eigen::VectorXd sqrt_weights = Eigen::Map<const Eigen::VectorXd>(weights.data(), weights.size()).array().sqrt();
    for (int i = 0; i < n; ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, i); it; ++it) {
            it.valueRef() *= sqrt_weights[i];
        }
        yvec[i] *= sqrt_weights[i];
        Eigen::VectorXd sqrt_weights = Eigen::Map<const Eigen::VectorXd>(weights.data(), weights.size()).array().sqrt();

    }
    // Build constraint mask: intercept unconstrained
    std::vector<bool> constrained(num_knots + 1, true);
    constrained[0] = false;

    // Solve
    auto box_start = high_resolution_clock::now();
    Eigen::VectorXd coeffs = box_lsq_ldlt(X, yvec, constrained);
    auto box_end = high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] box_lsq_ldlt sparse duration: " << std::chrono::duration<double>(box_end - box_start).count() << " seconds\n";

    // Predict and clamp
    Eigen::VectorXd yfit = X * coeffs;
    std::vector<double> result(yfit.size());
    std::transform(yfit.data(), yfit.data() + yfit.size(), result.begin(), [&](double v) {
        return clamp(v, min_val, max_val);
    });
    for (size_t i = 1; i < result.size(); ++i)
        result[i] = std::max(result[i], result[i-1]);

    auto end = high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] sparse fit_xy total duration: " << std::chrono::duration<double>(end - start).count() << " seconds\n";

    return result;
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
                    if (x[j] > 1e-12) {
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


InferPEP::InferPEP(bool use_ispline)
{
    if (use_ispline) {
        regressor_ptr_ = std::make_unique<IsplineRegression>();
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
