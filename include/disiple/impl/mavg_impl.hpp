#pragma once

#include <disiple/impl/fir_impl.hpp>
#include <disiple/impl/maybe_static.hpp>
#include <Eigen/Core>

namespace disiple {

    template <typename Scalar, int Length, int Stages>
    struct MAvgCoeffs
    : private MaybeStatic<Length, void>  // inheritance enables empty base-class optimization
    , private MaybeStatic<Stages, char>
    {
        explicit MAvgCoeffs(int n = Length, int s = Stages)
        : MaybeStatic<Length, void>(n)
        , MaybeStatic<Stages, char>(s)
        {}

        int length() const { return MaybeStatic<Length, void>::get(); }
        int stages() const { return MaybeStatic<Stages, char>::get(); }
    };

    template <int A, int B>
    struct multiply_extents
    {
        enum { value = A != Eigen::Dynamic && B != Eigen::Dynamic ? A*B : Eigen::Dynamic };
    };

    template <typename Scalar, int Length, int Channels, int Stages>
    struct MAvgState : FIRState<Scalar, Length, multiply_extents<Channels, Stages>::value>
    {
        enum { Rows = multiply_extents<Channels, Stages>::value };

        using Base = FIRState<Scalar, Length, Rows>;

        MAvgState() : num_(0), correct_num_(0)
        {
            if (Length != Eigen::Dynamic && Rows != Eigen::Dynamic)
                initialize();
        }

        void setup(const MAvgCoeffs<Scalar, Length, Stages>& coeffs, int nchans)
        {
            const Eigen::DenseIndex rws = coeffs.stages() * nchans;
            const Eigen::DenseIndex len = coeffs.length();
            if (Base::num_chans() != rws ||
                Base::length()    != len)
            {
                Base::buf_.resize(rws, len);
                sum_.resize(nchans, coeffs.stages());
                correct_sum_.resize(nchans, coeffs.stages());
                initialize();
            }
        }

        void initialize()
        {
            Base::initialize();
            sum_.setZero(); correct_sum_.setZero();
            num_ = 0;       correct_num_ = static_cast<int>(buf_.cols());
        }

        template <typename X>
        void apply(const MAvgCoeffs<Scalar, Length, Stages>& coeffs,
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

        using Base::advance;
        using Base::buf_;
        using Base::pos_;
        Eigen::Array<Scalar, Channels, Stages> sum_, correct_sum_;
        int num_, correct_num_;
    };


    // cumulative moving average

    template <typename Scalar>
    struct CumMAvgCoeffs
    {
        CumMAvgCoeffs() {}
        explicit CumMAvgCoeffs(int n) {}

        static int    length()  { return 1; }
    };

    template <typename Scalar, int Channels>
    struct CumMAvgState
    {
        CumMAvgState() : num_(0)
        {
            if (Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const CumMAvgCoeffs<Scalar>& coeffs, int nchans)
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
        void apply(const CumMAvgCoeffs<Scalar>&,
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
    };

}
