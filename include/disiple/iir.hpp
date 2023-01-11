#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <IIRImplementation N>
    struct Implementation { static constexpr IIRImplementation implementation = N; };

    template <typename... Options>
    using IIRParameters = Parameters<
        List<Options...>,
        OptionalValue<IIRImplementation, Implementation, DF2T>,
        OptionalValue<int, Stages, Eigen::Dynamic>,
        OptionalValue<int, Channels, 1>
    >;

    /// Infinite impulse response (IIR) digital filter
    template <typename Scalar, typename... Options>
    struct IIR : public FilterBase<Scalar, IIRParameters<Options...>::channels,
            IIRImpl<Scalar, IIRParameters<Options...>::stages, IIRParameters<Options...>::channels, IIRParameters<Options...>::implementation>,
            SecondOrderSections<Scalar, IIRParameters<Options...>::stages>
        >
    {
        using P      = IIRParameters<Options...>;
        using State  = IIRImpl<Scalar, P::stages, P::channels, P::implementation>;
        using Coeffs = SecondOrderSections<Scalar, P::stages>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        IIR() {}
        IIR(const IIRDesign &d) : Base(d) {}

        int num_stages() const { return Base::coeffs().num_stages(); }
        using Base::coeffs;
    };

}
