#pragma once

#include <disiple/impl/fir_impl.hpp>
#include <Eigen/Core>

namespace disiple {

    template <typename Scalar, int Length>
    struct DelayCoeffs
    {
        explicit DelayCoeffs(int n = Length)
        { assert(n == Length); }

        static int length() { return Length; }
    };

    template <typename Scalar>
    struct DelayCoeffs<Scalar, Eigen::Dynamic>
    {
        explicit DelayCoeffs(int n = 0)
        : length_(n) {}

        int length() const { return length_; }

        int length_;
    };

    template <typename Scalar, int Length, int Channels>
    struct DelayState : FIRState<Scalar, Length, Channels>
    {
        using Base = FIRState<Scalar, Length, Channels>;

        DelayState() : Base() {}

        using Base::initialize;
        using Base::advance;

        void setup(const DelayCoeffs<Scalar, Length>& coeffs, int nchans)
        {
            Base::setup(coeffs, nchans);
            prev_.resize(nchans);
        }

        template <typename X>
        void apply(const DelayCoeffs<Scalar, Length>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            advance();
            prev_ = buf_.col(pos_);
            buf_.col(pos_) = xi;
            xi = prev_;
        }

        template <typename X>
        void apply(const DelayCoeffs<Scalar, Length>& coeffs,
                   const Eigen::ArrayBase<X>& xi, DryRun)
        {
            advance();
            buf_.col(pos_) = xi;
        }

        Eigen::Array<Scalar, Channels, 1> prev_;
        using Base::buf_;
        using Base::pos_;
    };

}
