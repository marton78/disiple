#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/poly_fir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Polynomial FIR filter.
    /// A polynomial filter's weights are given by
    ///     $$ w_t = \sum_i a_i t^j $$
    /// where the $a_i$ are rational functions of the window length T:
    ///     $$ a_i = \sum_j p_ij n^j / \sum_k q_ik n^k $$
    /// with $t \in [1, T]$, $i \in [0, I)$, $j \in [0, J)$ and $k \in [0, K)$.
    template <typename Element, int I, int J, int K, int Length = Eigen::Dynamic>
    struct PolynomialFIR : public FilterBase<Element,
        PolyFIRState<typename ElementTraits<Element>::Scalar, Length, ElementTraits<Element>::Channels, I>,
        PolyFIRCoeffs<typename ElementTraits<Element>::Scalar, Length, I, J, K>
    >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = PolyFIRState<Scalar, Length, Channels, I>;
        using Coeffs = PolyFIRCoeffs<Scalar, Length, I, J, K>;
        using Base   = FilterBase<Element, State, Coeffs>;

        PolynomialFIR() {}

        /// Initialize the filter.
        template <typename P, typename Q>
        PolynomialFIR(const Eigen::ArrayBase<P>& p,   ///< Numerator of the weights, must be I x J
                      const Eigen::ArrayBase<Q>& q,   ///< Denominator of the weights, must be I x K
                      int window_length               ///< Length of the window. Must be equal to @Length if not dynamic.
        ) : Base(p, q, window_length)
        {
            static_assert(I != Eigen::Dynamic && J != Eigen::Dynamic && K != Eigen::Dynamic,
                          "All polynomials must be fixed-size!");
        }
    };
    
}
