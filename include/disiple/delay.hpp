#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/delay_impl.hpp>

namespace disiple {

    template <typename Scalar, int Length = Eigen::Dynamic, int Channels = 1>
    struct Delay : public FilterBase<Scalar, Channels,
        DelayState<Scalar, Length, Channels>,
        DelayCoeffs<Scalar, Length>
    >
    {
        using State  = DelayState<Scalar, Length, Channels>;
        using Coeffs = DelayCoeffs<Scalar, Length>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        Delay() {}
        explicit Delay(int length) : Base(length) {}
    };

}
