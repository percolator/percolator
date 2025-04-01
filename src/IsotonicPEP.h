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
    std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
                               double min_val, double max_val) const override {;};

private:
};


class IsplineRegression : public IsotonicRegression {
public:
    IsplineRegression(int degree = 2, int num_knots = 10)
        : degree_(degree), num_knots_(num_knots) {}

    std::vector<double> fit_y(const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override {
        int n = y.size();
        std::vector<double> x(n);
        std::iota(x.begin(), x.end(), 0);
        return fit_xy(x, y);
    }
        
    std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override {
        assert(x.size() == y.size());
        int n = x.size();
                
        create_knots(x);
                
        int num_basis = num_knots_ + degree_ - 1;
        Eigen::MatrixXd B(n, num_basis);
    
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < num_basis; ++j)
                B(i, j) = ispline_basis(j, degree_, x[i], knots_);
    
        Eigen::VectorXd yvec = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());
        Eigen::VectorXd coeffs = nnls(B, yvec);
    
        Eigen::VectorXd yfit = B * coeffs;
        return std::vector<double>(yfit.data(), yfit.data() + yfit.size());
    }
                
private:
    double mspline_basis(int i, int k, double x, const std::vector<double>& knots) const;
    double ispline_basis(int i, int k, double x, const std::vector<double>& knots, int steps = 50) const;
    void create_knots(const std::vector<double>& x) const;
    Eigen::VectorXd nnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const;

    int degree_;
    int num_knots_;
    mutable std::vector<double> knots_;
};

void IsplineRegression::create_knots(const std::vector<double>& x) const {
    knots_.resize(num_knots_ + 2 * degree_);
    for (size_t i = 0; i < knots_.size(); ++i) {
        if (i < degree_)
            knots_[i] = x.front();
        else if (i >= num_knots_ + degree_)
            knots_[i] = x.back();
        else {
            double alpha = double(i - degree_ + 1) / (num_knots_ + 1);
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

    for (int j = 1; j < steps; ++j) {
        double t = a + j * h;
        sum += mspline_basis(i, k, t, knots);
    }
    sum += 0.5 * mspline_basis(i, k, b, knots);
    return sum * h;
}

Eigen::VectorXd IsplineRegression::nnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const {
    Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
    Eigen::VectorXd w = A.transpose() * (b - A * x);

    std::vector<bool> passive_set(A.cols(), false);

    for (int iter = 0; iter < 500; ++iter) {
        int idx = -1;
        double max_w = 0;
        for (int i = 0; i < w.size(); ++i) {
            if (!passive_set[i] && w[i] > max_w) {
                max_w = w[i];
                idx = i;
            }
        }
        if (idx == -1) break;

        passive_set[idx] = true;
        std::vector<int> active;
        for (int i = 0; i < passive_set.size(); ++i)
            if (passive_set[i]) active.push_back(i);

        Eigen::MatrixXd A_active(A.rows(), active.size());
        for (int i = 0; i < active.size(); ++i)
            A_active.col(i) = A.col(active[i]);

        Eigen::VectorXd z = A_active.colPivHouseholderQr().solve(b);

        for (int i = 0; i < active.size(); ++i)
            x[active[i]] = std::max(0.0, z[i]);

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
            } else {
                regressor_ptr_ = std::make_unique<PavaRegression>();
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
            for (auto& dp : decoy_rate) {
                if (dp > 1. - epsilon)
                    dp = 1. - epsilon;
                double pep = dp / (1 - dp);
                if (pep > 1.)
                    pep = 1.;
                pep_iso.push_back(pep);
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