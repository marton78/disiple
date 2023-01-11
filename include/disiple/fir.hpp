#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/fir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Finite impulse response (FIR) digital filter
    template <typename Element, int Length = Eigen::Dynamic>
    struct FIR : public FilterBase<Element,
            FIRState<typename ElementTraits<Element>::Scalar, Length, ElementTraits<Element>::Channels>,
            FIRCoeffs<typename ElementTraits<Element>::Scalar, Length>
        >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = FIRState<Scalar, Length, Channels>;
        using Coeffs = FIRCoeffs<Scalar, Length>;
        using Base   = FilterBase<Element, State, Coeffs>;

        FIR() {}

        FIR(FIRDesign d) : Base(std::move(d)) {}

        template <typename A>
        FIR(const Eigen::ArrayBase<A>& a) : Base(a) {}

        int length() const { return Base::coeffs().length(); }
        using Base::coeffs;
    };

}
