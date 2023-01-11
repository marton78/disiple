#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/mavg_impl.hpp>

namespace disiple {

    template <typename Element, int Length = Eigen::Dynamic, int Stages = 1>
    struct moving_average : public filter_base<Element,
        mavg_state<typename element_traits<Element>::Scalar, Length, element_traits<Element>::Channels, Stages>,
        mavg_coeffs<typename element_traits<Element>::Scalar, Length, Stages>
    >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = mavg_state<Scalar, Length, Channels, Stages>;
        using Coeffs = mavg_coeffs<Scalar, Length, Stages>;
        using Base   = filter_base<Element, State, Coeffs>;

        moving_average() {}
        explicit moving_average(int length) : Base(length) {}
        moving_average(int length, int stages) : Base(length, stages) {}
    };


    template <typename Element>
    struct cum_moving_average : public filter_base<Element,
        cmavg_state<typename element_traits<Element>::Scalar, element_traits<Element>::Channels>,
        cmavg_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = cmavg_state<Scalar, Channels>;
        using Coeffs = cmavg_coeffs<Scalar>;
        using Base   = filter_base<Element, State, Coeffs>;

        cum_moving_average() {}
        explicit cum_moving_average(int length) : Base(length) {}
    };

}
