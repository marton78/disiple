#pragma once

#include <Eigen/Core>
#include <disiple/filter_design.hpp>

namespace disiple {

    template <typename Scalar, int Length>
    struct fir_coeffs
    {
        fir_coeffs() {
            if (Length != Eigen::Dynamic) {
                b_.setZero();
                b_[0] = 1;
            }
        }

        explicit fir_coeffs(fir_design d) : b_(d.coeffs.cast<Scalar>()) {}

        template <typename T>
        explicit fir_coeffs(T&& t) : b_(std::forward<T>(t)) {}

        int length() const { return static_cast<int>(b_.size()); }

        Eigen::Matrix<Scalar, Length, 1>  b_;
    };

    template <typename Scalar, int Length, int Channels>
    struct fir_state
    {
        fir_state() : pos_(0)
        {
            if (Channels != Eigen::Dynamic && Length != Eigen::Dynamic)
                initialize();
        }

        int num_chans() const { return static_cast<int>(buf_.rows()); }
        int length()    const { return static_cast<int>(buf_.cols()); }

        template <typename Coeffs>
        void setup(const Coeffs& coeffs, int nchans)
        {
            const int len = coeffs.length();
            if (num_chans() != nchans || length() != coeffs.length())
            {
                buf_.resize(nchans, len);
                initialize();
            }
        }

        void initialize()
        {
            buf_.fill(0); pos_ = 0;
        }

        void advance()
        {
            if (--pos_<0) pos_ += static_cast<int>(buf_.cols());
        }

        template <typename X>
        void apply(const fir_coeffs<Scalar, Length>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            advance();
            buf_.col(pos_) = xi;
            const int n0 = length() - pos_;
            const int n1 = pos_;
            xi = (buf_.rightCols(n0) * coeffs.b_.head(n0)).array();
            if (n1)
                xi += (buf_.leftCols(n1)  * coeffs.b_.tail(n1)).array();
        }

        template <typename X>
        void apply(const fir_coeffs<Scalar, Length>& coeffs,
                   const Eigen::ArrayBase<X>& xi, dry_run_t)
        {
            advance();
            buf_.col(pos_) = xi;
        }

        Eigen::Matrix<Scalar, Channels, Length> buf_;
        int pos_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

}
