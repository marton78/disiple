#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Infinite impulse response (IIR) digital filter
    template <typename Element, int Stages = Eigen::Dynamic, IIRImplementation Type = DF2T>
    struct iir : public filter_base<Element,
            iir_impl<typename element_traits<Element>::Scalar, Stages, element_traits<Element>::Channels, Type>,
            second_order_sections<typename element_traits<Element>::Scalar, Stages>
        >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = iir_impl<Scalar, Stages, Channels, Type>;
        using Coeffs = second_order_sections<Scalar, Stages>;
        using Base   = filter_base<Element, State, Coeffs>;

        iir() {}
        iir(const iir_design &d) : Base(d) {}

        int num_stages() const { return Base::coeffs().num_stages(); }
        using Base::coeffs;
    };

}
