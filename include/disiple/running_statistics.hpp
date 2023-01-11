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
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = running_minmax_state<Element, std::less>;
        using coeffs_type = running_minmax_coeffs<Scalar>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        explicit running_min(int len) : base_type(len) {}
    };

    template <typename Element>
    struct running_max : public filter_base<Element,
        running_minmax_state<Element, std::greater>,
        running_minmax_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = running_minmax_state<Element, std::greater>;
        using coeffs_type = running_minmax_coeffs<Scalar>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        explicit running_max(int len) : base_type(len) {}
    };

    template <typename Element>
    struct running_range : public filter_base<Element,
        running_range_state<Element>,
        running_minmax_coeffs<typename element_traits<Element>::Scalar>
    >
    {
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = running_range_state<Element>;
        using coeffs_type = running_minmax_coeffs<Scalar>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        explicit running_range(int len) : base_type(len) {}
    };

}
