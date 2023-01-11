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
        using Scalar = typename element_traits<Element>::Scalar;
        using State  = fir_state<Scalar, Length, Channels>;
        using Coeffs = fir_coeffs<Scalar, Length>;
        using Base   = filter_base<Element, State, Coeffs>;

        fir() {}

        fir(fir_design d) : Base(std::move(d)) {}

        template <typename A>
        fir(const Eigen::ArrayBase<A>& a) : Base(a) {}

        int length() const { return Base::coeffs().length(); }
        using Base::coeffs;
    };

}
