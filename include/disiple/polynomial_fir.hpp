#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/poly_fir_impl.hpp>
#include <disiple/filter_design.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    /// Polynomial FIR filter.
    /// A polynomial filter's weights are given by
    ///     $$ w_t = \sum_i a_i t^j $$
    /// where the $a_i$ are rational functions of the window length T:
    ///     $$ a_i = \sum_j p_ij n^j / \sum_k q_ik n^k $$
    /// with $t \in [1, T]$, $i \in [0, I)$, $j \in [0, J)$ and $k \in [0, K)$.
    template <typename Scalar, int I, int J, int K, typename... Options>
    struct PolynomialFIR : public FilterBase<Scalar, PolynomialFIR<Scalar, I, J, K, Options...>>,
                           public Parameters<
                                List<Options...>,
                                OptionalValue<int, Length, Eigen::Dynamic>,
                                OptionalValue<int, Channels, 1>
                            >
    {
    public:
        using State  = PolyFIRState<Scalar, PolynomialFIR::LengthValue::Static, PolynomialFIR::ChannelsValue::Static, I>;
        using Coeffs = PolyFIRCoeffs<Scalar, PolynomialFIR::LengthValue::Static, I, J, K>;

        PolynomialFIR() {}

        /// Initialize the filter.
        template <typename P, typename Q>
        PolynomialFIR(const Eigen::ArrayBase<P>& p,   ///< Numerator of the weights, must be I x J
                      const Eigen::ArrayBase<Q>& q,   ///< Denominator of the weights, must be I x K
                      int window_length               ///< Length of the window. Must be equal to @Length if not dynamic.
        ) : coeffs_(p, q, window_length)
        {
            static_assert(I != Eigen::Dynamic && J != Eigen::Dynamic && K != Eigen::Dynamic,
                          "All polynomials must be fixed-size!");
        }

    private:
        friend FilterBase<Scalar, PolynomialFIR<Scalar, I, J, K, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };
    
}
