#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/mavg_impl.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <typename... Options>
    using MovingAverageParameters = Parameters<
        List<Options...>,
        OptionalValue<int, Length, Eigen::Dynamic>,
        OptionalValue<int, Stages, 1>,
        OptionalValue<int, Channels, 1>
    >;

    template <typename Scalar, typename... Options>
    struct MovingAverage : public FilterBase<Scalar, MovingAverageParameters<Options...>::channels,
        MAvgState <Scalar, MovingAverageParameters<Options...>::length, MovingAverageParameters<Options...>::channels, MovingAverageParameters<Options...>::stages>,
        MAvgCoeffs<Scalar, MovingAverageParameters<Options...>::length, MovingAverageParameters<Options...>::stages>
    >
    {
        using P      = MovingAverageParameters<Options...>;
        using State  = MAvgState <Scalar, P::length, P::channels, P::stages>;
        using Coeffs = MAvgCoeffs<Scalar, P::length, P::stages>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        MovingAverage() {}
        explicit MovingAverage(int length) : Base(length) {}
        MovingAverage(int length, int stages) : Base(length, stages) {}
    };


    template <typename Scalar, int Channels = 1>
    struct CumMovingAverage : public FilterBase<Scalar, Channels,
        CumMAvgState<Scalar, Channels>,
        CumMAvgCoeffs<Scalar>
    >
    {
        using State  = CumMAvgState<Scalar, Channels>;
        using Coeffs = CumMAvgCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        CumMovingAverage() {}
        explicit CumMovingAverage(int length) : Base(length) {}
    };

}
