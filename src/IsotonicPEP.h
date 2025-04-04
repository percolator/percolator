#ifndef ISOTONIC_REGRESSION_H_
#define ISOTONIC_REGRESSION_H_

#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>    // for std::iota
#include <cassert>
#include <iostream>
#include <Eigen/Dense>

template <typename T>
const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

class IsotonicRegression {
public:
    IsotonicRegression() = default;
    virtual ~IsotonicRegression() = default;

    virtual std::vector<double> fit_y(const std::vector<double>& y, 
                                      double min_val = std::numeric_limits<double>::lowest(),
                                      double max_val = std::numeric_limits<double>::max()) const = 0;

    virtual std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
                                       double min_val = std::numeric_limits<double>::lowest(),
                                       double max_val = std::numeric_limits<double>::max()) const = 0;
};

class PavaRegression : public IsotonicRegression {
protected:
    //--------------------------------------------------------
    // A simple block structure for merging
    //--------------------------------------------------------
    struct Block {
        double sum   = 0.0;
        int    count = 0;
        double avg   = 0.0;
    };

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
    virtual std::vector<double> pavaNonDecreasingRanged(
        const std::vector<double>& values,
        const double min_value = std::numeric_limits<double>::min(),
        const double max_value = std::numeric_limits<double>::max()
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
            double val = max( min(b.avg, max_value), min_value);
            for (int c = 0; c < b.count; ++c) {
                result.push_back(val);
            }
        }
        return result;
    }

public:
    std::vector<double> fit_y(const std::vector<double>& y, double min_val, double max_val) const override
    { return pavaNonDecreasingRanged(y, min_val, max_val); }
    std::vector<double> fit_xy(const std::vector<double>&/* x */, const std::vector<double>& y,
                               double min_val, double max_val) const override { return fit_y(y, min_val, max_val); };

};


class IsplineRegression : public IsotonicRegression {
public:
    IsplineRegression(int degree = 2)
        : degree_(degree) {}
    
    std::vector<double> fit_y(const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override {
        int n = y.size();
        std::vector<double> x(n);
        std::iota(x.begin(), x.end(), 0);
        return fit_xy(x, y, min_val, max_val);
    }
        
    std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
        double min_val = std::numeric_limits<double>::lowest(),
        double max_val = std::numeric_limits<double>::max()) const override {
    
        assert(x.size() == y.size());
        int n = x.size();
        if (n == 0) return {};

        // Normalize x to [0,1]
        std::vector<double> x_norm(n);
        double xmin = *std::min_element(x.begin(), x.end());
        double xmax = *std::max_element(x.begin(), x.end());
        double xscale = xmax - xmin;
        for (int i = 0; i < n; ++i)
            x_norm[i] = (x[i] - xmin) / xscale;

        // Determine number of knots and construct them
        int num_knots = std::min(static_cast<int>(std::sqrt(n)), 500);
        std::vector<double> knots(num_knots + 1);
        for (int i = 0; i <= num_knots; ++i)
            knots[i] = static_cast<double>(i) / num_knots;

        // Build cubic I-spline basis
        Eigen::MatrixXd B(n, num_knots);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < num_knots; ++j)
                B(i, j) = cubic_ispline(x_norm[i], knots[j], knots[j+1]);

        // Add intercept column
        Eigen::MatrixXd X(n, num_knots + 1);
        X.col(0) = Eigen::VectorXd::Ones(n);
        X.block(0, 1, n, num_knots) = B;

        // Build constraint mask: intercept (c[0]) is unconstrained, rest must be ≥ 0
        std::vector<bool> constrained(num_knots + 1, true);
        constrained[0] = false;

        // Solve using LDLT with constraint projection
        Eigen::VectorXd yvec = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());
        Eigen::VectorXd coeffs = box_lsq_ldlt(X, yvec, constrained);

        // Predict and clamp
        Eigen::VectorXd yfit = X * coeffs;
        std::vector<double> result(yfit.size());
        std::transform(yfit.data(), yfit.data() + yfit.size(), result.begin(), [&](double v) {
            return clamp(v, min_val, max_val);
        });

        return result;
    }
    /*
    std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override {
        assert(x.size() == y.size());
        int n = x.size();

        int num_knots = min(static_cast<int>(sqrt(n)),500);
                
        create_knots(x, num_knots);
                
        int num_basis = num_knots + degree_ - 1;
        Eigen::MatrixXd B(n, num_basis);
    
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < num_basis; ++j)
                B(i, j) = cubic_ispline(x[i], knots_[j], knots_[j+1]);
//                B(i, j) = ispline_basis(j, degree_, x[i], knots_);
    
        Eigen::VectorXd yvec = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());
        Eigen::VectorXd coeffs = lsqnonneg(B, yvec);
    
        Eigen::VectorXd yfit = B * coeffs;
        std::vector<double> result(yfit.size());
        // auto ret = std::vector<double>(yfit.data(), yfit.data() + yfit.size());
        std::transform(yfit.data(), yfit.data() + yfit.size(), result.begin(), [&](double v) {
            return clamp(v, min_val, max_val);
        });
        return result;
    }
    */            
private:
    double mspline_basis(int i, int k, double x, const std::vector<double>& knots) const;
    double ispline_basis(int i, int k, double x, const std::vector<double>& knots, int steps = 50) const;
    void create_knots(const std::vector<double>& x, const int& num_knots) const;
    double cubic_ispline(double x, double left, double right) const;
    // Eigen::VectorXd nnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const;
    Eigen::VectorXd lsqnonneg(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const;
    Eigen::VectorXd box_lsq_ldlt(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const std::vector<bool>& constrained,
        double lambda = 0.0, int max_iter = 100, double tol = 1e-12) const;
    size_t degree_;
    mutable std::vector<double> knots_;
};

void IsplineRegression::create_knots(const std::vector<double>& x, const int& num_knots) const {
    knots_.resize(num_knots + 2 * degree_);
    for (size_t i = 0; i < knots_.size(); ++i) {
        if (i < degree_)
            knots_[i] = x.front();
        else if (i >= num_knots + degree_)
            knots_[i] = x.back();
        else {
            double alpha = double(i - degree_ + 1) / (num_knots + 1);
            knots_[i] = x.front() + alpha * (x.back() - x.front());
        }
    }
}

double IsplineRegression::mspline_basis(int i, int k, double x, const std::vector<double>& t) const {
    if (k == 0) {
        return (t[i] <= x && x < t[i+1]) ? 1.0 / (t[i+1] - t[i]) : 0.0;
    }
    double denom = t[i + k + 1] - t[i];
    double left = 0.0, right = 0.0;

    if (t[i + k] > t[i])
        left = (x - t[i]) * mspline_basis(i, k - 1, x, t) / (t[i + k] - t[i]);
    if (t[i + k + 1] > t[i + 1])
        right = (t[i + k + 1] - x) * mspline_basis(i + 1, k - 1, x, t) / (t[i + k + 1] - t[i + 1]);

    return k * (left + right) / denom;
}

double IsplineRegression::ispline_basis(int i, int k, double x, const std::vector<double>& knots, int steps) const {
    double a = knots.front();
    double b = x;
    double h = (b - a) / steps;
    double sum = 0.5 * mspline_basis(i, k, a, knots);

    #pragma omp parallel for reduction(+:sum)
    for (int j = 1; j < steps; ++j) {
        double t = a + j * h;
        sum += mspline_basis(i, k, t, knots);
    }
    sum += 0.5 * mspline_basis(i, k, b, knots);
    return sum * h;
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

Eigen::VectorXd IsplineRegression::lsqnonneg(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const {
    double tol = 1e-10;

    using namespace Eigen;
    const int m = A.rows();
    const int n = A.cols();

    VectorXd x = VectorXd::Zero(n);
    VectorXd w = A.transpose() * (b - A * x);

    std::vector<bool> P(n, false);
    std::vector<int> passive_set;

    for (int iter = 0; iter < 3 * n; ++iter) {
        // Step 1: find index with max w
        int t = -1;
        double max_w = 0;
        for (int i = 0; i < n; ++i) {
            if (!P[i] && w[i] > max_w) {
                max_w = w[i];
                t = i;
            }
        }
        if (t == -1) break;

        P[t] = true;
        passive_set.push_back(t);

        VectorXd z = x;
        bool inner_loop = true;
        while (inner_loop) {
            // Solve least squares on passive set
            MatrixXd A_p(m, passive_set.size());
            for (size_t i = 0; i < passive_set.size(); ++i)
                A_p.col(i) = A.col(passive_set[i]);

            VectorXd b_p = A_p.colPivHouseholderQr().solve(b);

            z.setZero();
            for (size_t i = 0; i < passive_set.size(); ++i)
                z[passive_set[i]] = b_p[i];

            if ((z.array() >= 0).all()) {
                x = z;
                inner_loop = false;
            } else {
                double alpha = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < passive_set.size(); ++i) {
                    int idx = passive_set[i];
                    if (z[idx] < 0)
                        alpha = std::min(alpha, x[idx] / (x[idx] - z[idx]));
                }

                x = x + alpha * (z - x);

                std::vector<int> new_passive;
                for (int idx : passive_set) {
                    if (x[idx] > tol)
                        new_passive.push_back(idx);
                    else
                        P[idx] = false;
                }
                passive_set = std::move(new_passive);
            }
        }
        w = A.transpose() * (b - A * x);
    }

    return x;
}


class InferPEP {
    public:
        InferPEP(bool use_ispline)
        {
            if (use_ispline) {
                regressor_ptr_ = std::make_unique<IsplineRegression>();
                if (VERB > 1) {
                    cerr << "Performing isotonic regression using I-Splines" << endl;
                }                
            } else {
                regressor_ptr_ = std::make_unique<PavaRegression>();
                if (VERB > 1) {
                    cerr << "Performing isotonic regression using PAVA" << endl;
                }
            }
        }
    
        std::vector<double> q_to_pep(const std::vector<double>& q_values) {
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
    
        std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
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
    
        std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {}) {
            double epsilon = 1e-20;
            auto is_dec = is_decoy;
            is_dec.insert(is_dec.begin(), 0.5);
    
            std::vector<double> decoy_rate;
            if (!scores.empty()) {
                auto sc = scores;
                sc.insert(sc.begin(), sc[0]);
                decoy_rate = regressor_ptr_->fit_xy(sc, is_dec);
            } else {
                decoy_rate = regressor_ptr_->fit_y(is_dec);
            }
            decoy_rate.erase(decoy_rate.begin());
    
            std::vector<double> pep_iso;
            auto i_d = is_decoy.begin(); double tar = 0.; double dec = 0.; double errors = 0.;
            for (auto& dp : decoy_rate) {
                if (dp > 1. - epsilon)
                    dp = 1. - epsilon;
                double pep = dp / (1 - dp);
                if (pep > 1.)
                    pep = 1.;
                pep_iso.push_back(pep);
                if (*i_d>0.5) {dec+=1.;} else {tar+=1.; errors += pep;} i_d++;
                cerr << dp << "\t" << pep << "\t" << tar << "\t" << dec << "\t" << errors << "\t" << (dec+0.5)/(tar) << "\t" << errors/(tar) << endl;
            }
            return pep_iso;
        }

        double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
            double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
            return interp_pep;
        }
        
    private:
        std::unique_ptr<IsotonicRegression> regressor_ptr_;
        std::vector<double> qs;
        std::vector<double> pep_iso;
    };

    #endif /* ISOTONICPEP_H_ */