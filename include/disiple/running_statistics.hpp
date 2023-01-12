#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/running_stats_impl.hpp>
#include <disiple/named_params.hpp>
#include <deque>

namespace disiple {

    /// Filter to calculate running minimum.
    /// Algorithm from: David Lemire, “Streaming Maximum-Minimum Filter Using No
    /// More than Three Comparisons per Element”,
    /// Nordic Journal of Computing,  Vol. 13, 2006

    template <typename Scalar, typename... Options>
    class RunningMin : public FilterBase<Scalar, RunningMin<Scalar, Options...>>,
                       public Parameters<
                                List<Options...>,
                                OptionalValue<int, Channels, 1>
                            >
    {
    public:
        using State  = RunningMinMaxState<Scalar, RunningMin::ChannelsValue::Static, std::less>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        friend FilterBase<Scalar, RunningMin<Scalar, Options...>>;

        explicit RunningMin(int len) : coeffs_(len) {}

    private:
        friend FilterBase<Scalar, RunningMin<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

    template <typename Scalar, typename... Options>
    class RunningMax : public FilterBase<Scalar, RunningMax<Scalar, Options...>>,
                       public Parameters<
                                List<Options...>,
                                OptionalValue<int, Channels, 1>
                            >
    {
    public:
        using State  = RunningMinMaxState<Scalar, RunningMax::ChannelsValue::Static, std::greater>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        friend FilterBase<Scalar, RunningMax<Scalar, Options...>>;

        explicit RunningMax(int len) : coeffs_(len) {}

    private:
        friend FilterBase<Scalar, RunningMax<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

    template <typename Scalar, typename... Options>
    class RunningRange : public FilterBase<Scalar, RunningRange<Scalar, Options...>>,
                         public Parameters<
                                List<Options...>,
                                OptionalValue<int, Channels, 1>
                            >
    {
    public:
        using State  = RunningRangeState<Scalar, RunningRange::ChannelsValue::Static>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        friend FilterBase<Scalar, RunningRange<Scalar, Options...>>;

        explicit RunningRange(int len) : coeffs_(len) {}

    private:
        friend FilterBase<Scalar, RunningRange<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

}
