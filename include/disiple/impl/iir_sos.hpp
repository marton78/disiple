#pragma once

#include <disiple/filter_design.hpp>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

namespace disiple {

    template <typename Scalar, int Stages = Eigen::Dynamic>
    class SecondOrderSections
    {
    public:
        using Coeffs = Eigen::Array<Scalar, 4, Stages>;
        using Column = typename Coeffs::ConstColXpr;

        SecondOrderSections() : scaling_(1) { coeffs_.setZero(); }
        SecondOrderSections(const IIRDesign& l);

        std::complex<Scalar> response(Scalar f) const;

        template <typename F, typename Z>
        void response(const Eigen::ArrayBase<F>& f, Eigen::ArrayBase<Z>& z) const;

        int num_stages() const { return int(coeffs_.cols()); }

        // The coefficients of the filter are:
        // b = [1 b_1 b_2] * scaling;
        // a = [1 a_1 a_2]; with a_i = -m_i;
        Scalar b1(int i) const { return coeffs_(0,i); }
        Scalar b2(int i) const { return coeffs_(1,i); }
        Scalar m1(int i) const { return coeffs_(2,i); }
        Scalar m2(int i) const { return coeffs_(3,i); }
        Scalar scaling() const { return scaling_; }

        Column coeffs(int i) const { return coeffs_.col(i); }

    protected:
        Coeffs  coeffs_;
        Scalar    scaling_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    struct wrong_static_size : std::runtime_error
    {
        wrong_static_size(std::string msg) : std::runtime_error(std::move(msg)) {}
    };

    namespace internal {

        using namespace Eigen;

        template <typename Scalar>
        void biquad_from_pz_pair(Ref<Array<Scalar, 4, 1>, Aligned> coeffs,
                                 Complex pole1, Complex zero1,
                                 Complex pole2, Complex zero2);

        template <typename Scalar>
        void response(Ref<const Array<Scalar, 4, Dynamic>, Aligned, Stride<4, 1>> coeffs,
                      Scalar                                                      scaling,
                      Ref<const Array<Scalar, Dynamic, 1>>                        f,
                      Ref<Array<std::complex<Scalar>, Dynamic, 1>>                z
                      );

    }

    template <typename Scalar, int Stages>
    SecondOrderSections<Scalar, Stages>::SecondOrderSections(const IIRDesign& l) : scaling_(1)
    {
        const int n = (l.num_poles()+1)/2;

        if (Stages != Eigen::Dynamic && Stages != n)
            throw wrong_static_size(std::string("Static size is wrong: n=") + std::to_string(n));

        coeffs_.resize(4, n);

        for (int i=0; i<n; ++i)
            internal::biquad_from_pz_pair<Scalar>(coeffs_.col(i),
                                                  l[i].poles.first,  l[i].zeros.first,
                                                  l[i].poles.second, l[i].zeros.second);

        scaling_ = Scalar(l.normal_gain()) / std::abs(response(Scalar(l.normal_w())));
    }

    template <typename Scalar, int Stages>
    std::complex<Scalar> SecondOrderSections<Scalar, Stages>::response(Scalar f) const
    {
        using Complex = std::complex<Scalar>;

        const Complex czn1(std::cos(-f), std::sin(-f));
        const Complex czn2(std::cos(-Scalar(2)*f), std::sin(-Scalar(2)*f));

        Complex z(scaling_);

        for (int i=0; i<num_stages(); ++i)
        {
            auto c = coeffs_.col(i);
            z *= (Complex(1) + c[0] * czn1 + c[1] * czn2)
               / (Complex(1) - c[2] * czn1 - c[3] * czn2);
        }

        return z;
    }

    template <typename Scalar, int Stages>
    template <typename F, typename Z>
    void SecondOrderSections<Scalar, Stages>::response(const Eigen::ArrayBase<F>& f, Eigen::ArrayBase<Z>& z) const
    {
        internal::response<Scalar>(coeffs_, scaling_, f.derived(), z.derived());
    }

}
