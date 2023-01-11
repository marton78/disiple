#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/running_stats_impl.hpp>
#include <deque>

namespace disiple {

    /// Filter to calculate running minimum.
    /// Algorithm from: David Lemire, “Streaming Maximum-Minimum Filter Using No
    /// More than Three Comparisons per Element”,
    /// Nordic Journal of Computing,  Vol. 13, 2006

    template <typename Element>
    struct RunningMin : public FilterBase<Element,
        RunningMinMaxState<Element, std::less>,
        RunningMinMaxCoeffs<typename ElementTraits<Element>::Scalar>
    >
    {
        using Scalar = typename ElementTraits<Element>::Scalar;
        using State  = RunningMinMaxState<Element, std::less>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Element, State, Coeffs>;

        explicit RunningMin(int len) : Base(len) {}
    };

    template <typename Element>
    struct RunningMax : public FilterBase<Element,
        RunningMinMaxState<Element, std::greater>,
        RunningMinMaxCoeffs<typename ElementTraits<Element>::Scalar>
    >
    {
        using Scalar      = typename ElementTraits<Element>::Scalar;
        using State  = RunningMinMaxState<Element, std::greater>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Element, State, Coeffs>;

        explicit RunningMax(int len) : Base(len) {}
    };

    template <typename Element>
    struct RunningRange : public FilterBase<Element,
        RunningRangeState<Element>,
        RunningMinMaxCoeffs<typename ElementTraits<Element>::Scalar>
    >
    {
        using Scalar      = typename ElementTraits<Element>::Scalar;
        using State  = RunningRangeState<Element>;
        using Coeffs = RunningMinMaxCoeffs<Scalar>;
        using Base   = FilterBase<Element, State, Coeffs>;

        explicit RunningRange(int len) : Base(len) {}
    };

}
