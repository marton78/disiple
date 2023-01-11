#pragma once

#include <disiple/impl/fir_impl.hpp>
#include <Eigen/Core>

namespace disiple {

    template <typename Scalar, int Length>
    struct delay_coeffs
    {
        explicit delay_coeffs(int n = Length)
        { assert(n == Length); }

        static int length() { return Length; }
    };

    template <typename Scalar>
    struct delay_coeffs<Scalar, Eigen::Dynamic>
    {
        explicit delay_coeffs(int n = 0)
        : length_(n) {}

        int length() const { return length_; }

        int length_;
    };

    template <typename Scalar, int Length, int Channels>
    struct delay_state : fir_state<Scalar, Length, Channels>
    {
        using Base = fir_state<Scalar, Length, Channels>;

        delay_state() : Base() {}

        using Base::initialize;
        using Base::advance;

        void setup(const delay_coeffs<Scalar, Length>& coeffs, int nchans)
        {
            Base::setup(coeffs, nchans);
            prev_.resize(nchans);
        }

        template <typename X>
        void apply(const delay_coeffs<Scalar, Length>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            advance();
            prev_ = buf_.col(pos_);
            buf_.col(pos_) = xi;
            xi = prev_;
        }

        template <typename X>
        void apply(const delay_coeffs<Scalar, Length>& coeffs,
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
