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
    struct running_min : public filter_base<Element,
        running_minmax_state<Element, std::less>,
        running_minmax_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = running_minmax_state<Element, std::less>;
        using Coeffs = running_minmax_coeffs<Scalar>;
        using Base   = filter_base<Element, State, Coeffs>;

        explicit running_min(int len) : Base(len) {}
    };

    template <typename Element>
    struct running_max : public filter_base<Element,
        running_minmax_state<Element, std::greater>,
        running_minmax_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        using Scalar      = typename element_traits<Element>::Scalar;
        using State  = running_minmax_state<Element, std::greater>;
        using Coeffs = running_minmax_coeffs<Scalar>;
        using Base   = filter_base<Element, State, Coeffs>;

        explicit running_max(int len) : Base(len) {}
    };

    template <typename Element>
    struct running_range : public filter_base<Element,
        running_range_state<Element>,
        running_minmax_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        using Scalar      = typename element_traits<Element>::Scalar;
        using State  = running_range_state<Element>;
        using Coeffs = running_minmax_coeffs<Scalar>;
        using Base   = filter_base<Element, State, Coeffs>;

        explicit running_range(int len) : Base(len) {}
    };

}
