#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Infinite impulse response (IIR) digital filter
    template <typename Element, int Stages = Eigen::Dynamic, implementation_type Type = DF2T>
    struct iir : public filter_base<Element,
            iir_impl<typename element_traits<Element>::Scalar, Stages, element_traits<Element>::Channels, Type>,
            second_order_sections<typename element_traits<Element>::Scalar, Stages>
        >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = iir_impl<Scalar, Stages, Channels, Type>;
        using coeffs_type = second_order_sections<Scalar, Stages>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        iir() {}
        iir(const iir_design &d) : base_type(d) {}

        int num_stages() const { return base_type::coeffs().num_stages(); }
        using base_type::coeffs;
    };

}
