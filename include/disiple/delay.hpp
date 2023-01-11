#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/delay_impl.hpp>

namespace disiple {

    template <typename Element, int Length = Eigen::Dynamic>
    struct delay : public filter_base<Element,
        delay_state<typename element_traits<Element>::Scalar, Length, element_traits<Element>::Channels>,
        delay_coeffs<typename element_traits<Element>::Scalar, Length>
    >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = delay_state<Scalar, Length, Channels>;
        using coeffs_type = delay_coeffs<Scalar, Length>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        delay() {}
        explicit delay(int length) : base_type(length) {}
    };

}
