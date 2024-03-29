#pragma once

#include <Eigen/Core>
#include <disiple/impl/mavg_impl.hpp>

namespace disiple {

    // polynomial fir filter

    template <typename Scalar, int Length, int Stages, int J, int K>
    struct PolyFIRCoeffs : MAvgCoeffs<Scalar, Length, Stages>
    {
        using Base = MAvgCoeffs<Scalar, Length, Stages>;

        PolyFIRCoeffs() {}

        template <typename P, typename Q>
        explicit PolyFIRCoeffs(const Eigen::ArrayBase<P>& p, const Eigen::ArrayBase<Q>& q, int wlen) :
            Base(wlen, p.rows()), p_(p), q_(q)
        {}

        template <typename A>
        void update(Scalar n, Eigen::ArrayBase<A>& a) const
        {
            Eigen::Array<Scalar, Stages, 1> sum_p = p_.col(0), sum_q = q_.col(0);

            enum { M = J<K ? J : K };
            Scalar nn = Scalar(1);
            for (int i = 1; i < M; ++i) { nn *= n; sum_p += nn * p_.col(i); sum_q += nn * q_.col(i); }
            for (int j = M; j < J; ++j) { nn *= n; sum_p += nn * p_.col(j); }
            for (int k = M; k < K; ++k) { nn *= n; sum_q += nn * q_.col(k); }
            a = (sum_q == 0).select(Scalar(0), sum_p / sum_q);
        }

        Eigen::Array<Scalar, Stages, J>  p_;
        Eigen::Array<Scalar, Stages, K>  q_;
    };

    template <typename Scalar, int Length, int Channels, int Stages>
    struct PolyFIRState : FIRState<Scalar, Length, Channels>
    {
        using Base = FIRState<Scalar, Length, Channels>;

        PolyFIRState() : num_(0)
        {
            if (Length != Eigen::Dynamic && Channels != Eigen::Dynamic)
                initialize();
        }

        template <int J, int K>
        void setup(const PolyFIRCoeffs<Scalar, Length, Stages, J, K>& coeffs, int nchans)
        {
            const int len = coeffs.length();
            if (Base::num_chans() != nchans ||
                Base::length()    != len)
            {
                Base::buf_.resize(nchans, len);
                sum_.resize(nchans, coeffs.stages());
                a_.resize(coeffs.stages());
                initialize();
            }
        }

        void initialize()
        {
            Base::initialize();
            sum_.fill(0); num_ = 0; a_.fill(0);
        }

        template <int J, int K, typename X>
        void apply(const PolyFIRCoeffs<Scalar, Length, Stages, J, K>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            advance();
            auto x_n = buf_.col(pos_).array();

            if (num_ == buf_.cols())
            {
                // remove oldest sample
                Scalar n = Scalar(num_);
                sum_.col(0) -= x_n;
                for (int i = 1; i < coeffs.stages(); ++i)
                    sum_.col(i) -= (x_n *= n);
            }
            else
            {
                // re-calculate coefficients
                Scalar n = Scalar(++num_);
                coeffs.update(n, a_);
            }

            // store new sample
            x_n = xi;

            // update state and calculate output
            xi = a_[0] * (sum_.col(0) += x_n);

            for (int i = 1; i < coeffs.stages(); ++i)
                xi += a_[i] * (sum_.col(i) += sum_.col(i-1));
        }

        using Base::advance;
        using Base::buf_;
        using Base::pos_;
        Eigen::Array<Scalar, Channels, Stages>      sum_;
        Eigen::Array<Scalar, Stages, 1>             a_;
        int                                         num_;
    };

}
