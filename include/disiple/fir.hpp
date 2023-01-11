#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/fir_impl.hpp>
#include <disiple/filter_design.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <typename... Options>
    using FIRParameters = Parameters<
        List<Options...>,
        OptionalValue<int, Length, Eigen::Dynamic>,
        OptionalValue<int, Channels, 1>
    >;

    /// Finite impulse response (FIR) digital filter
    template <typename Scalar, typename... Options>
    struct FIR : public FilterBase<Scalar, FIRParameters<Options...>::channels,
            FIRState <Scalar, FIRParameters<Options...>::length, FIRParameters<Options...>::channels>,
            FIRCoeffs<Scalar, FIRParameters<Options...>::length>
        >
    {
        using P      = FIRParameters<Options...>;
        using State  = FIRState<Scalar, P::length, P::channels>;
        using Coeffs = FIRCoeffs<Scalar, P::length>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        FIR() {}

        FIR(FIRDesign d) : Base(std::move(d)) {}

        template <typename A>
        FIR(const Eigen::ArrayBase<A>& a) : Base(a) {}

        int length() const { return Base::coeffs().length(); }
        using Base::coeffs;
    };

}
