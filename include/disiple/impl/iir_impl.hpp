#pragma once

#include <disiple/impl/iir_sos.hpp>
#include <Eigen/Core>

namespace disiple {

    enum implementation_type { DF1, DF2, DF1T, DF2T };

    template <typename Scalar, int Stages, int Channels, implementation_type Type>
    struct iir_impl;


    // Direct Form I

    template <typename Scalar, int Stages, int Channels>
    struct iir_impl<Scalar, Stages, Channels, DF1>
    {
        iir_impl()
        {
            if (Stages != Eigen::Dynamic && Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const second_order_sections<Scalar, Stages>& sos, int nchans)
        {
            const int nstages = sos.num_stages();
            if (history_.rows() != nchans || history_.cols() != 4 * nstages)
            {
                history_.resize(nchans, 4 * nstages);
                initialize();
            }
        }

        void initialize()
        {
            history_.fill(0);
        }

        template <typename X>
        void initialize(const second_order_sections<Scalar, Stages>& sos,
                        const Eigen::ArrayBase<X>& x_ss)
        {
            Scalar r = (Scalar(1) + sos.b1(0) + sos.b2(0))
                     / (Scalar(1) - sos.m1(0) - sos.m2(0));
            history_.col(0) = history_.col(1) = x_ss;
            history_.col(2) = history_.col(3) = history_.col(0) * r;

            for (int j=1, k=4; k<history_.cols(); ++j, k+=4)
            {
                auto h = history_.template block<Channels, 4>(0, k, history_.rows(), 4);

                r = (Scalar(1) + sos.b1(j) + sos.b2(j))
                  / (Scalar(1) - sos.m1(j) - sos.m2(j));
                h.col(0) = h.col(1) = history_.col(k-2);  // x1 = x2 = y1;
                h.col(2) = h.col(3) = h.col(0) * r;       // y1 = y2 = x1 * r;
            }
        }

        template <typename X>
        void apply(const second_order_sections<Scalar, Stages>& sos,
                   Eigen::ArrayBase<X>& xi)
        {
            for (int j=0, k=0; k<history_.cols(); ++j, k+=4)
            {
                // direct form 1:
                //     y = x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;

                auto h = history_.template block<Channels, 4>(0, k, history_.rows(), 4);

                x0_ = xi;
                xi += (h.matrix() * sos.coeffs(j).matrix()).array();

                h.col(1) = h.col(0);  // x2 = x1, x1 = xi
                h.col(0) = x0_;
                h.col(3) = h.col(2);  // y2 = y1, y1 = yi
                h.col(2) = xi;
            }
            xi *= sos.scaling();
        }

        enum { S4 = Stages != Eigen::Dynamic ? 4 * Stages : Eigen::Dynamic };

        Eigen::Array<Scalar, Channels, 1>   x0_;
        Eigen::Array<Scalar, Channels, S4>  history_; // x1 x2 y1 y2

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };




    // Direct Form II

    template <typename Scalar, int Stages, int Channels>
    struct iir_impl<Scalar, Stages, Channels, DF2>
    {
        iir_impl()
        {
            if (Stages != Eigen::Dynamic && Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const second_order_sections<Scalar, Stages>& sos, int nchans)
        {
            const int nstages = sos.num_stages();
            if (w1_.rows() != nchans || w1_.cols() != nstages)
            {
                w1_.resize(nchans, nstages); w2_.resize(nchans, nstages);
                initialize();
            }
        }

        void initialize()
        {
            w1_.fill(0); w2_.fill(0);
        }

        template <typename X>
        void initialize(const second_order_sections<Scalar, Stages>& sos,
                        const Eigen::ArrayBase<X>& x_ss)
        {
            Scalar r = Scalar(1) / (Scalar(1) - sos.m1(0) - sos.m2(0));
            w1_.col(0) = w2_.col(0) = x_ss * r;

            for (int j=1; j<w1_.cols(); ++j)
            {
                r = (Scalar(1) + sos.b1(j-1) + sos.b2(j-1))
                  / (Scalar(1) - sos.m1(j)   - sos.m2(j)  );
                w1_.col(j) = w2_.col(j) = w1_.col(0) * r;
            }
        }

        template <typename X>
        void apply(const second_order_sections<Scalar, Stages>& sos,
                   Eigen::ArrayBase<X>& xi)
        {
            for (int j=0; j<w1_.cols(); ++j)
            {
                w0_ = xi  + sos.m1(j) * w1_.col(j) + sos.m2(j) * w2_.col(j);
                xi  = w0_ + sos.b1(j) * w1_.col(j) + sos.b2(j) * w2_.col(j);

                w2_.col(j) = w1_.col(j);
                w1_.col(j) = w0_;
            }
            xi *= sos.scaling();
        }

        Eigen::Array<Scalar, Channels, 1>      w0_;
        Eigen::Array<Scalar, Channels, Stages> w1_, w2_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };





    // Direct Form II Transposed

    template <typename Scalar, int Stages, int Channels>
    struct iir_impl<Scalar, Stages, Channels, DF2T>
    {
        iir_impl()
        {
            if (Stages != Eigen::Dynamic && Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const second_order_sections<Scalar, Stages>& sos, int nchans)
        {
            const int nstages = sos.num_stages();
            if (d1_.rows() != nchans || d1_.cols() != nstages)
            {
                d1_.resize(nchans, nstages); d2_.resize(nchans, nstages);
                initialize();
            }
        }

        void initialize()
        {
            d1_.fill(0); d2_.fill(0);
        }

        template <typename X>
        void initialize(const second_order_sections<Scalar, Stages>& sos,
                        const Eigen::ArrayBase<X>& x_ss)
        {
            Scalar r = (Scalar(1) + sos.b1(0) + sos.b2(0))
                     / (Scalar(1) - sos.m1(0) - sos.m2(0));

            d1_.col(0) = x_ss;
            d2_.col(0) = (sos.b2(0) + sos.m2(0) * r) * d1_.col(0);

            for (int j=1; j<d1_.cols(); ++j)
            {
                d1_.col(j) = d1_.col(j-1) * r;
                d1_.col(j-1) *= (r - Scalar(1));

                r = (Scalar(1) + sos.b1(j) + sos.b2(j))
                  / (Scalar(1) - sos.m1(j) - sos.m2(j));

                d2_.col(j) = (sos.b2(j) + sos.m2(j) * r) * d1_.col(j);
            }
            d1_.col(d1_.cols()-1) *= (r - Scalar(1));
        }

        template <typename X>
        void apply(const second_order_sections<Scalar, Stages>& sos,
                   Eigen::ArrayBase<X>& xi)
        {
            for (int j=0; j<d1_.cols(); ++j)
            {
                x0_ = xi;
                xi += d1_.col(j);
                d1_.col(j) = sos.b1(j) * x0_ + sos.m1(j) * xi + d2_.col(j);
                d2_.col(j) = sos.b2(j) * x0_ + sos.m2(j) * xi;
            }
            xi *= sos.scaling();
        }

        Eigen::Array<Scalar, Channels, 1>      x0_;
        Eigen::Array<Scalar, Channels, Stages> d1_, d2_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

}
