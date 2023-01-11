#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/running_stats_impl.hpp>
#include <disiple/named_params.hpp>
#include <deque>

namespace disiple {

    template <typename... Options>
    using RunningStatsParameters = Parameters<
        List<Options...>,
        OptionalValue<int, Channels, 1>
    >;

    /// Filter to calculate running minimum.
    /// Algorithm from: David Lemire, “Streaming Maximum-Minimum Filter Using No
    /// More than Three Comparisons per Element”,
    /// Nordic Journal of Computing,  Vol. 13, 2006

    template <typename Scalar, typename... Options>
    struct RunningMin : public FilterBase<Scalar, RunningStatsParameters<Options...>::channels,
        RunningMinMaxState<Scalar, RunningStatsParameters<Options...>::channels, std::less>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using P      = RunningStatsParameters<Options...>;
        using State  = RunningMinMaxState<Scalar, P::channels, std::less>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        explicit RunningMin(int len) : Base(len) {}
    };

    template <typename Scalar, typename... Options>
    struct RunningMax : public FilterBase<Scalar, RunningStatsParameters<Options...>::channels,
        RunningMinMaxState<Scalar, RunningStatsParameters<Options...>::channels, std::greater>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using P      = RunningStatsParameters<Options...>;
        using State  = RunningMinMaxState<Scalar, P::channels, std::greater>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        explicit RunningMax(int len) : Base(len) {}
    };

    template <typename Scalar, typename... Options>
    struct RunningRange : public FilterBase<Scalar, RunningStatsParameters<Options...>::channels,
        RunningRangeState<Scalar, RunningStatsParameters<Options...>::channels>,
        RunningMinMaxCoeffs<Scalar>
    >
    {
        using P      = RunningStatsParameters<Options...>;
        using State  = RunningRangeState<Scalar, P::channels>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        explicit RunningRange(int len) : Base(len) {}
    };

}
