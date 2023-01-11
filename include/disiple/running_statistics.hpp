#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/running_stats_impl.hpp>
#include <deque>

namespace disiple {

    /// Filter to calculate running minimum.
    /// Algorithm from: David Lemire, “Streaming Maximum-Minimum Filter Using No
    /// More than Three Comparisons per Element”,
    /// Nordic Journal of Computing,  Vol. 13, 2006

    template <typename Scalar, int Channels = 1>
    struct RunningMin : public FilterBase<Scalar, Channels,
        RunningMinMaxState<Scalar, Channels, std::less>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using State  = RunningMinMaxState<Scalar, Channels, std::less>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        explicit RunningMin(int len) : Base(len) {}
    };

    template <typename Scalar, int Channels = 1>
    struct RunningMax : public FilterBase<Scalar, Channels,
        RunningMinMaxState<Scalar, Channels, std::greater>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using State  = RunningMinMaxState<Scalar, Channels, std::greater>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        explicit RunningMax(int len) : Base(len) {}
    };

    template <typename Scalar, int Channels = 1>
    struct RunningRange : public FilterBase<Scalar, Channels,
        RunningRangeState<Scalar, Channels>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using State  = RunningRangeState<Scalar, Channels>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, Channels, State, Coeffs>;

        explicit RunningRange(int len) : Base(len) {}
    };

}
