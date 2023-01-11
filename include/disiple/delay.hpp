#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/delay_impl.hpp>
#include <disiple/named_params.hpp>


namespace disiple {

    template <typename... Options>
    using DelayParameters = Parameters<
        List<Options...>,
        OptionalValue<int, Length, Eigen::Dynamic>,
        OptionalValue<int, Channels, 1>
    >;

    template <typename Scalar, typename... Options>
    struct Delay : public FilterBase<Scalar, DelayParameters<Options...>::channels,
        DelayState <Scalar, DelayParameters<Options...>::length, DelayParameters<Options...>::channels>,
        DelayCoeffs<Scalar, DelayParameters<Options...>::length>
    >
    {
        using P      = DelayParameters<Options...>;
        using State  = DelayState<Scalar, P::length, P::channels>;
        using Coeffs = DelayCoeffs<Scalar, P::length>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        Delay() {}
        explicit Delay(int length) : Base(length) {}
    };

}
