#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/fir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Finite impulse response (FIR) digital filter
    template <typename Scalar, int Length = Eigen::Dynamic, int Channels = 1>
    struct FIR : public FilterBase<Scalar, Channels,
            FIRState<Scalar, Length, Channels>,
            FIRCoeffs<Scalar, Length>
        >
    {
        using State  = FIRState<Scalar, Length, Channels>;
        using Coeffs = FIRCoeffs<Scalar, Length>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        FIR() {}

        FIR(FIRDesign d) : Base(std::move(d)) {}

        template <typename A>
        FIR(const Eigen::ArrayBase<A>& a) : Base(a) {}

        int length() const { return Base::coeffs().length(); }
        using Base::coeffs;
    };

}
