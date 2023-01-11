#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/mavg_impl.hpp>

namespace disiple {

    template <typename Scalar, int Length = Eigen::Dynamic, int Stages = 1, int Channels = 1>
    struct MovingAverage : public FilterBase<Scalar, Channels,
        MAvgState<Scalar, Length, Channels, Stages>,
        MAvgCoeffs<Scalar, Length, Stages>
    >
    {
        using State  = MAvgState<Scalar, Length, Channels, Stages>;
        using Coeffs = MAvgCoeffs<Scalar, Length, Stages>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

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
