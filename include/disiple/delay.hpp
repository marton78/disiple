#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/delay_impl.hpp>

namespace disiple {

    template <typename Element, int Length = Eigen::Dynamic>
    struct Delay : public FilterBase<Element,
        DelayState<typename ElementTraits<Element>::Scalar, Length, ElementTraits<Element>::Channels>,
        DelayCoeffs<typename ElementTraits<Element>::Scalar, Length>
    >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = DelayState<Scalar, Length, Channels>;
        using Coeffs = DelayCoeffs<Scalar, Length>;
        using Base   = FilterBase<Element, State, Coeffs>;

        Delay() {}
        explicit Delay(int length) : Base(length) {}
    };

}
