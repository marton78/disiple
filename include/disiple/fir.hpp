#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/fir_impl.hpp>
#include <disiple/filter_design.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    /// Finite impulse response (FIR) digital filter
    template <typename Scalar, typename... Options>
    class FIR : public FilterBase<Scalar, FIR<Scalar, Options...>>,
                public Parameters<
                        List<Options...>,
                        OptionalValue<int, Length, Eigen::Dynamic>,
                        OptionalValue<int, Channels, 1>
                    >
    {
    public:
        using State = FIRState<Scalar, FIR::length, FIR::channels>;
        using Coeffs = FIRCoeffs<Scalar, FIR::length>;

        FIR() {}

        FIR(FIRDesign d) : coeffs_(std::move(d)) {}

        template <typename A>
        FIR(const Eigen::ArrayBase<A>& a) : coeffs_(a) {}

        int length() const { return coeffs_.length(); }

    private:
        friend FilterBase<Scalar, FIR<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

}
