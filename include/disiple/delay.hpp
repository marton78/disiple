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
        typedef typename element_traits<Element>::Scalar        Scalar;
        typedef delay_state<Scalar, Length, Channels>           state_type;
        typedef delay_coeffs<Scalar, Length>                    coeffs_type;
        typedef filter_base<Element, state_type, coeffs_type>   base_type;

        delay() {}
        explicit delay(int length) : base_type(length) {}
    };

}
