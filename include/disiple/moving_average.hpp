#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/mavg_impl.hpp>

namespace disiple {

    template <typename Element, int Length = Eigen::Dynamic, int Stages = 1>
    struct MovingAverage : public FilterBase<Element,
        MAvgState<typename ElementTraits<Element>::Scalar, Length, ElementTraits<Element>::Channels, Stages>,
        MAvgCoeffs<typename ElementTraits<Element>::Scalar, Length, Stages>
    >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = MAvgState<Scalar, Length, Channels, Stages>;
        using Coeffs = MAvgCoeffs<Scalar, Length, Stages>;
        using Base   = FilterBase<Element, State, Coeffs>;

        MovingAverage() {}
        explicit MovingAverage(int length) : Base(length) {}
        MovingAverage(int length, int stages) : Base(length, stages) {}
    };


    template <typename Element>
    struct CumMovingAverage : public FilterBase<Element,
        CumMAvgState<typename ElementTraits<Element>::Scalar, ElementTraits<Element>::Channels>,
        CumMAvgCoeffs<typename ElementTraits<Element>::Scalar>
    >
    {
        enum { Channels = ElementTraits<Element>::Channels };
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = CumMAvgState<Scalar, Channels>;
        using Coeffs = CumMAvgCoeffs<Scalar>;
        using Base   = FilterBase<Element, State, Coeffs>;

        CumMovingAverage() {}
        explicit CumMovingAverage(int length) : Base(length) {}
    };

}
