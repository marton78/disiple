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
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = delay_state<Scalar, Length, Channels>;
        using Coeffs = delay_coeffs<Scalar, Length>;
        using Base   = filter_base<Element, State, Coeffs>;

        delay() {}
        explicit delay(int length) : Base(length) {}
    };

}
