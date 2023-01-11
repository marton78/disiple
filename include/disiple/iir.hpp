#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Infinite impulse response (IIR) digital filter
    template <typename Scalar, int Stages = Eigen::Dynamic, IIRImplementation Type = DF2T, int Channels = 1>
    struct IIR : public FilterBase<Scalar, Channels,
            IIRImpl<Scalar, Stages, Channels, Type>,
            SecondOrderSections<Scalar, Stages>
        >
    {
        using State  = IIRImpl<Scalar, Stages, Channels, Type>;
        using Coeffs = SecondOrderSections<Scalar, Stages>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        IIR() {}
        IIR(const IIRDesign &d) : Base(d) {}

        int num_stages() const { return Base::coeffs().num_stages(); }
        using Base::coeffs;
    };

}
