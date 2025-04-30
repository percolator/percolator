#ifndef ISOTONIC_REGRESSION_H_
#define ISOTONIC_REGRESSION_H_

#include <chrono>
#include <iostream> // for std::cerr, std::endl
#include <memory>   // for std::unique_ptr, std::make_unique
#include <algorithm> // for std::min, std::max
#include <Eigen/Dense>
#include <Eigen/Sparse> 

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

    virtual std::vector<double> pavaNonDecreasingRanged(
        const std::vector<double>& values,
        const double min_value = std::numeric_limits<double>::min(),
        const double max_value = std::numeric_limits<double>::max()
    ) const;

public:
    std::vector<double> fit_y(const std::vector<double>& y, double min_val, double max_val) const override
    { return pavaNonDecreasingRanged(y, min_val, max_val); }
    std::vector<double> fit_xy(const std::vector<double>&/* x */, const std::vector<double>& y,
                               double min_val, double max_val) const override { return fit_y(y, min_val, max_val); };

};

struct BinnedData {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> weights;
};


class IsplineRegression : public IsotonicRegression {
public:
    IsplineRegression(int degree = 2)
        : degree_(degree) {}
    
    std::vector<double> fit_y(const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override;
        
    std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
        double min_val = std::numeric_limits<double>::lowest(),
        double max_val = std::numeric_limits<double>::max()) const override;

    double cubic_ispline(double x, double left, double right) const;

protected:
    std::vector<double> normalize_vector(const std::vector<double>& x, double xmin, double xscale) const;
    std::vector<double> compute_knots(const std::vector<double>& x_norm, const std::vector<double>& y, int num_knots) const;
    Eigen::SparseMatrix<double, Eigen::RowMajor> build_design_matrix(const std::vector<double>& x_norm,
        const std::vector<double>& knots, int n, int num_knots) const;
    void fill_intercept(Eigen::SparseMatrix<double, Eigen::RowMajor>& X) const;
    void apply_weights(Eigen::SparseMatrix<double, Eigen::RowMajor>& X,
        Eigen::VectorXd& yvec, const Eigen::VectorXd& sqrt_weights_vec) const;
    Eigen::MatrixXd build_dense_matrix(const std::vector<double>& x_full_norm,
        const std::vector<double>& knots, int n_full, int num_knots) const;
    BinnedData bin_data_weighted(const std::vector<double>& x, const std::vector<double>& y, int max_bins = 2500) const;
    // Eigen::VectorXd nnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) const;
    Eigen::VectorXd box_lsq_ldlt(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const std::vector<bool>& constrained,
        double lambda = 0.0, int max_iter = 100, double tol = 1e-12) const;
    size_t degree_;
    mutable std::vector<double> knots_;
};




class InferPEP {
    public:
        InferPEP(bool use_ispline);
    
        std::vector<double> q_to_pep(const std::vector<double>& q_values);
        std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores);
        std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {});

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