#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Infinite impulse response (IIR) digital filter
    template <typename Element, int Stages = Eigen::Dynamic, IIRImplementation Type = DF2T>
    struct IIR : public FilterBase<Element,
            IIRImpl<typename ElementTraits<Element>::Scalar, Stages, ElementTraits<Element>::Channels, Type>,
            SecondOrderSections<typename ElementTraits<Element>::Scalar, Stages>
        >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = IIRImpl<Scalar, Stages, Channels, Type>;
        using Coeffs = SecondOrderSections<Scalar, Stages>;
        using Base   = FilterBase<Element, State, Coeffs>;

        IIR() {}
        IIR(const IIRDesign &d) : Base(d) {}

        int num_stages() const { return Base::coeffs().num_stages(); }
        using Base::coeffs;
    };

}
