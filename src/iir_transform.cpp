#include <disiple/filter_design.hpp>

namespace disiple {

    static const complex_t infinity(std::numeric_limits<double>::infinity());

    /*
     * s-plane to z-plane transforms
     *
     * For pole filters, an analog prototype is created via placement of
     * poles and zeros in the s-plane. The analog prototype is either
     * a halfband low pass or a halfband low shelf. The poles, zeros,
     * and normalization parameters are transformed into the z-plane
     * using variants of the bilinear transformation.
     *
     */

    //------------------------------------------------------------------------------

    void lowpass::operator()(const iir_prototype& analog, iir_design& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;
        const double f = tan(pi_half * cutoff_);

        auto transform = [f] (complex_t c) -> complex_t
        {
            if (c == infinity)
                return complex_t(-1, 0);

            // frequency transform
            c = f * c;

            // bilinear low pass transform
            return (1. + c) / (1. - c);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            digital.add_conjugate_pair(transform(pair.poles.first), transform(pair.zeros.first));
        }

        if (numPoles & 1)
        {
            const pole_zero_pair& pair = analog[pairs];
            assert(pair.poles.first.imag()==0);
            assert(pair.zeros.first.imag()==0);
            complex_t p = transform(pair.poles.first);
            complex_t z = transform(pair.zeros.first);
            digital.add_single(p.real(), z.real());
        }

        digital.set_normal(analog.normal_w(), analog.normal_gain());
    }


    //------------------------------------------------------------------------------

    void highpass::operator()(const iir_prototype& analog, iir_design& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;
        const double f = 1. / tan(pi_half * cutoff_);

        auto transform = [f] (complex_t c) -> complex_t
        {
            if (c == infinity)
                return complex_t (1, 0);

            // frequency transform
            c = f * c;

            // bilinear high pass transform
            return - (1. + c) / (1. - c);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            digital.add_conjugate_pair(transform(pair.poles.first),
                                       transform(pair.zeros.first));
        }

        if (numPoles & 1)
        {
            const pole_zero_pair& pair = analog[pairs];
            assert(pair.poles.first.imag()==0);
            assert(pair.zeros.first.imag()==0);
            complex_t p = transform(pair.poles.first);
            complex_t z = transform(pair.zeros.first);
            digital.add_single(p.real(), z.real());
        }

        digital.set_normal(pi - analog.normal_w(), analog.normal_gain());
    }

    //------------------------------------------------------------------------------

    void bandpass::operator()(const iir_prototype& analog, iir_design& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;

        double wc2 = lo_ * pi_half;
        double wc  = hi_ * pi_half;

        // what is this crap?
        if (wc2 < 1e-8)
            wc2 = 1e-8;
        if (wc  > pi_half-1e-8)
            wc  = pi_half-1e-8;

        const double a = cos(wc + wc2) / cos(wc - wc2);
        const double b = 1.0 / tan(wc - wc2);
        const double ab_2 = 2.0 * a*b;
        const double t = b*b*(a*a-1.0);
        const double t1 = 4.0 * (t+1);
        const double t2 = 8.0 * (t-1);

        auto transform = [&] (complex_t c) -> complex_pair
        {
            if (c == infinity)
                return complex_pair(-1, 1);

            c = (1. + c) / (1. - c); // bilinear

            complex_t w = ab_2 * c + ab_2;
            complex_t v = std::sqrt(((t1*c + t2) * c) + t1);
            complex_t d = 2.0 * ((b-1)*c + (b+1));

            return complex_pair((w-v)/d, (w+v)/d);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            complex_pair p1 = transform(pair.poles.first);
            complex_pair z1 = transform(pair.zeros.first);

            //
            // Optimize out the calculations for conjugates for Release builds
            //
#if !defined(NDEBUG)
            complex_pair p2 = transform(pair.poles.second);
            assert (p2.first == std::conj(p1.first));
            assert (p2.second == std::conj(p1.second));
#endif

            digital.add_conjugate_pair(p1.first, z1.first);
            digital.add_conjugate_pair(p1.second, z1.second);
        }

        if (numPoles & 1)
        {
            const pole_zero_pair& pair = analog[pairs];
            assert(pair.poles.first.imag()==0);
            assert(pair.zeros.first.imag()==0);
            complex_pair poles = transform(pair.poles.first);
            complex_pair zeros = transform(pair.zeros.first);

            digital.add(poles, zeros);
        }

        double wn = analog.normal_w() * 0.5;
        digital.set_normal(2 * atan(sqrt(tan(wc + wn) * tan(wc2 + wn))),
                           analog.normal_gain());
    }

    //------------------------------------------------------------------------------

    void bandstop::operator()(const iir_prototype& analog, iir_design& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;

        double wc2 = lo_ * pi_half;
        double wc  = hi_ * pi_half;

        // this is crap
        if (wc2 < 1e-8)
            wc2 = 1e-8;
        if (wc  > pi-1e-8)
            wc  = pi-1e-8;

        const double a = cos(wc + wc2) / cos(wc - wc2);
        const double b = tan(wc - wc2);
        const double a2 = a * a;
        const double b2 = b * b;
        const double t1 = 4.0 * (b2 + a2 - 1.0);
        const double t2 = 8.0 * (b2 - a2 + 1.0);

        auto transform = [&] (complex_t c) -> complex_pair
        {
            if (c == infinity)
                c = -1;
            else
                c = (1. + c) / (1. - c); // bilinear

            complex_t w = a * (1.0-c);
            complex_t u = 0.5 * std::sqrt((t1 * c + t2) * c + t1);
            complex_t d = (b-1.0)*c + (b+1.0);

            return complex_pair((w+u)/d, (w-u)/d);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            complex_pair p  = transform(pair.poles.first);
            complex_pair z  = transform(pair.zeros.first);

            //
            // Optimize out the calculations for conjugates for Release builds
            //
//#ifdef NDEBUG
//            // trick to get the conjugate
//            if (z.second == z.first)
//                z.second = std::conj (z.first);
//
//#else
//            // Do the full calculation to verify correctness
//            //            complex_pair pc = transform (analog[i].poles.second);
//            complex_pair zc = transform (analog[i].zeros.second);
//
//            // get the conjugates into pc and zc
//            if (zc.first == z.first)
//                std::swap (zc.first, zc.second);
//
//            //            assert (pc.first  == std::conj (p.first));
//            //            assert (pc.second == std::conj (p.second));
//            //            assert (zc.first  == std::conj (z.first));
//            //            assert (zc.second == std::conj (z.second));
//
//#endif
            digital.add_conjugate_pair (p.first, z.first);
            digital.add_conjugate_pair (p.second, z.second);
        }

        if (numPoles & 1)
        {
            const pole_zero_pair& pair = analog[pairs];
            assert(pair.poles.first.imag()==0);
            assert(pair.zeros.first.imag()==0);

            complex_pair poles = transform(pair.poles.first);
            complex_pair zeros = transform(pair.zeros.first);

            digital.add(poles, zeros);
        }

        if (wc+wc2 < pi_half)
            digital.set_normal(pi, analog.normal_gain());
        else
            digital.set_normal(0, analog.normal_gain());
    }

}
