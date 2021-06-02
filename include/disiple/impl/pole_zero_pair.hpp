#pragma once

#include <complex>
#include <cassert>
#include <utility>

namespace disiple {

    static const double pi      = 3.1415926535897932384626433832795028841971;
    static const double two_pi  = pi * 2.0;
    static const double pi_half = pi * 0.5;
    typedef std::complex<double> complex_t;

    struct complex_pair : std::pair<complex_t, complex_t>
    {
        complex_pair() {}

        explicit complex_pair(const complex_t& c1)
        : std::pair<complex_t, complex_t>(c1, 0.)
        { assert(c1.imag() == 0); }

        complex_pair(const complex_t& c1, const complex_t& c2)
        : std::pair<complex_t, complex_t>(c1, c2)
        {}

        bool is_conjugate() const { return second == std::conj(first); }

        // Returns true if this is either a conjugate pair,
        // or a pair of reals where neither is zero.
        bool is_matched_pair() const
        {
            if (first.imag() != 0)
                return second == std::conj(first);
            else
                return second.imag() == 0 &&
                second.real() != 0 &&
                first.real() != 0;
        }
    };


    // A pair of pole/zeros. This fits in a biquad (but is missing the gain)
    struct pole_zero_pair
    {
        complex_pair poles;
        complex_pair zeros;

        pole_zero_pair() {}

        // single pole/zero
        pole_zero_pair(const complex_t& p, const complex_t& z)
        : poles(p), zeros(z)
        {}

        // pole/zero pair
        pole_zero_pair(const complex_t& pole1, const complex_t& zero1,
                       const complex_t& pole2, const complex_t& zero2)
        : poles(pole1, pole2) , zeros(zero1, zero2)
        {}

        bool is_single_pole() const
        {
            return poles.second == 0. && zeros.second == 0.;
        }

    };

}
