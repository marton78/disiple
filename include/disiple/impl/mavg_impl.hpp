#pragma once

#include <disiple/impl/fir_impl.hpp>
#include <Eigen/Core>

namespace disiple {

    template <typename Scalar, int Length, int Stages>
    struct mavg_coeffs
    {
        explicit mavg_coeffs(int n = Length, int s = Stages)
        { assert(n == Length && s == Stages); }

        static int    length()  { return Length; }
        static int    stages()  { return Stages; }
    };

    template <typename Scalar, int Stages>
    struct mavg_coeffs<Scalar, Eigen::Dynamic, Stages>
    {
        explicit mavg_coeffs(int n = 0, int s = Stages)
        : length_(n) { assert(s == Stages); }

        int           length()  const { return length_; }
        static int    stages()        { return Stages; }

        int length_;
    };

    template <typename Scalar, int Length>
    struct mavg_coeffs<Scalar, Length, Eigen::Dynamic>
    {
        explicit mavg_coeffs(int n = Length, int s = 0)
        : stages_(s) { assert(n == Length); }

        static int    length()        { return Length; }
        int           stages()  const { return stages_; }

        int stages_;
    };

    template <typename Scalar>
    struct mavg_coeffs<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    {
        explicit mavg_coeffs(int n = 0, int s = 0)
        : length_(n), stages_(s) {}

        int           length()  const { return length_; }
        int           stages()  const { return stages_; }

        int length_, stages_;
    };

    template <int A, int B>
    struct multiply_extents
    {
        enum { value = A != Eigen::Dynamic && B != Eigen::Dynamic ? A*B : Eigen::Dynamic };
    };

    template <typename Scalar, int Length, int Channels, int Stages>
    struct mavg_state : fir_state<Scalar, Length, multiply_extents<Channels, Stages>::value>
    {
        enum { Rows = multiply_extents<Channels, Stages>::value };

        typedef fir_state<Scalar, Length, Rows> base_type;

        mavg_state() : num_(0), correct_num_(0)
        {
            if (Length != Eigen::Dynamic && Rows != Eigen::Dynamic)
                initialize();
        }

        void setup(const mavg_coeffs<Scalar, Length, Stages>& coeffs, int nchans)
        {
            const Eigen::DenseIndex rws = coeffs.stages() * nchans;
            const Eigen::DenseIndex len = coeffs.length();
            if (base_type::num_chans() != rws ||
                base_type::length()    != len)
            {
                base_type::buf_.resize(rws, len);
                sum_.resize(nchans, coeffs.stages());
                correct_sum_.resize(nchans, coeffs.stages());
                initialize();
            }
        }

        void initialize()
        {
            base_type::initialize();
            sum_.setZero(); correct_sum_.setZero();
            num_ = 0;       correct_num_ = static_cast<int>(buf_.cols());
        }

        template <typename X>
        void apply(const mavg_coeffs<Scalar, Length, Stages>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            using namespace Eigen;

            // Workaround: apparently Eigen doesn't align Array<Scalar, 1, 1>
            enum { A = Rows != 1 ? Aligned : Unaligned };

            Map<Array<Scalar, Rows, 1>, A> vsum(sum_.data(), buf_.rows(), 1);

            advance();

            if (num_ == buf_.cols())
                vsum -= buf_.col(pos_).array();
            else
                ++num_;

            const auto nchans = sum_.rows();
            const Scalar num  = Scalar(num_);
            for (int i = 0; i < coeffs.stages(); ++i)
            {
                buf_.template block<Channels, 1>(i * nchans, pos_, nchans, 1) = xi;
                correct_sum_.col(i) += xi;
                sum_.col(i) += xi;
                xi = sum_.col(i) / num;
            }

            if (--correct_num_ == 0) {
                sum_ = correct_sum_;
                correct_sum_.setZero();
                correct_num_ = static_cast<int>(buf_.cols());
            }

        }

        using base_type::advance;
        using base_type::buf_;
        using base_type::pos_;
        Eigen::Array<Scalar, Channels, Stages> sum_, correct_sum_;
        int num_, correct_num_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };


    // cumulative moving average

    template <typename Scalar>
    struct cmavg_coeffs
    {
        cmavg_coeffs() {}
        explicit cmavg_coeffs(int n) {}

        static int    length()  { return 1; }
    };

    template <typename Scalar, int Channels>
    struct cmavg_state
    {
        cmavg_state() : num_(0)
        {
            if (Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const cmavg_coeffs<Scalar>& coeffs, int nchans)
        {
            if (avg_.size() != nchans)
            {
                avg_.resize(nchans);
                initialize();
            }
        }

        void initialize()
        {
            avg_.fill(0); num_ = 0;
        }

        template <typename X>
        void apply(const cmavg_coeffs<Scalar>&,
                   Eigen::ArrayBase<X>& xi)
        {
            ++num_;
            if (num_ == 1) {
                xi.fill(0);
            } else {
                avg_ += (xi - avg_) / Scalar(num_);
                xi = avg_;
            }
        }

        Eigen::Array<Scalar, Channels, 1>      avg_;
        int num_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

}
