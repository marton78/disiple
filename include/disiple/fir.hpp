#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/fir_impl.hpp>
#include <disiple/filter_design.hpp>

namespace disiple {

    /// Finite impulse response (FIR) digital filter
    template <typename Element, int Length = Eigen::Dynamic>
    struct fir : public filter_base<Element,
            fir_state<typename element_traits<Element>::Scalar, Length, element_traits<Element>::Channels>,
            fir_coeffs<typename element_traits<Element>::Scalar, Length>
        >
    {
        enum { Channels = element_traits<Element>::Channels };
        using Scalar      = typename element_traits<Element>::Scalar;
        using state_type  = fir_state<Scalar, Length, Channels>;
        using coeffs_type = fir_coeffs<Scalar, Length>;
        using base_type   = filter_base<Element, state_type, coeffs_type>;

        fir() {}

        fir(fir_design d) : base_type(std::move(d)) {}

        template <typename A>
        fir(const Eigen::ArrayBase<A>& a) : base_type(a) {}

        int length() const { return base_type::coeffs().length(); }
        using base_type::coeffs;
    };

}
