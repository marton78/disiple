#include <disiple/filter_design.hpp>

namespace disiple {

    static const Complex infinity(std::numeric_limits<double>::infinity());

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

    void Lowpass::operator()(const IIRPrototype& analog, IIRDesign& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;
        const double f = tan(pi_half * cutoff_);

        auto transform = [f] (Complex c) -> Complex
        {
            if (c == infinity)
                return Complex(-1, 0);

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
            Complex p = transform(pair.poles.first);
            Complex z = transform(pair.zeros.first);
            digital.add_single(p.real(), z.real());
        }

        digital.set_normal(analog.normal_w(), analog.normal_gain());
    }


    //------------------------------------------------------------------------------

    void Highpass::operator()(const IIRPrototype& analog, IIRDesign& digital) const
    {
        const int numPoles = analog.num_poles();
        const int pairs = numPoles / 2;
        const double f = 1. / tan(pi_half * cutoff_);

        auto transform = [f] (Complex c) -> Complex
        {
            if (c == infinity)
                return Complex (1, 0);

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
            Complex p = transform(pair.poles.first);
            Complex z = transform(pair.zeros.first);
            digital.add_single(p.real(), z.real());
        }

        digital.set_normal(pi - analog.normal_w(), analog.normal_gain());
    }

    //------------------------------------------------------------------------------

    void Bandpass::operator()(const IIRPrototype& analog, IIRDesign& digital) const
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

        auto transform = [&] (Complex c) -> ComplexPair
        {
            if (c == infinity)
                return ComplexPair(-1, 1);

            c = (1. + c) / (1. - c); // bilinear

            Complex w = ab_2 * c + ab_2;
            Complex v = std::sqrt(((t1*c + t2) * c) + t1);
            Complex d = 2.0 * ((b-1)*c + (b+1));

            return ComplexPair((w-v)/d, (w+v)/d);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            ComplexPair p1 = transform(pair.poles.first);
            ComplexPair z1 = transform(pair.zeros.first);

            //
            // Optimize out the calculations for conjugates for Release builds
            //
#if !defined(NDEBUG)
            ComplexPair p2 = transform(pair.poles.second);
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
            ComplexPair poles = transform(pair.poles.first);
            ComplexPair zeros = transform(pair.zeros.first);

            digital.add(poles, zeros);
        }

        double wn = analog.normal_w() * 0.5;
        digital.set_normal(2 * atan(sqrt(tan(wc + wn) * tan(wc2 + wn))),
                           analog.normal_gain());
    }

    //------------------------------------------------------------------------------

    void Bandstop::operator()(const IIRPrototype& analog, IIRDesign& digital) const
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

        auto transform = [&] (Complex c) -> ComplexPair
        {
            if (c == infinity)
                c = -1;
            else
                c = (1. + c) / (1. - c); // bilinear

            Complex w = a * (1.0-c);
            Complex u = 0.5 * std::sqrt((t1 * c + t2) * c + t1);
            Complex d = (b-1.0)*c + (b+1.0);

            return ComplexPair((w+u)/d, (w-u)/d);
        };

        for (int i = 0; i < pairs; ++i)
        {
            const pole_zero_pair& pair = analog[i];
            ComplexPair p  = transform(pair.poles.first);
            ComplexPair z  = transform(pair.zeros.first);

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
//            //            ComplexPair pc = transform (analog[i].poles.second);
//            ComplexPair zc = transform (analog[i].zeros.second);
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

            ComplexPair poles = transform(pair.poles.first);
            ComplexPair zeros = transform(pair.zeros.first);

            digital.add(poles, zeros);
        }

        if (wc+wc2 < pi_half)
            digital.set_normal(pi, analog.normal_gain());
        else
            digital.set_normal(0, analog.normal_gain());
    }

}
